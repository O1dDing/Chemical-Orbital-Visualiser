#include "cov/nbo_channels.hpp"
#include "cov/local_angular_projection.hpp"
#include "cov/orbital_chemistry.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <tuple>
namespace cov { namespace {
using V=std::array<double,3>;
V cross(V a,V b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
V unit(V a){double n=std::hypot(a[0],a[1],a[2]);if(n<1e-14)return {};for(auto& x:a)x/=n;return a;}
std::array<double,9> frame(const Atom& a,const Atom& b){
    const auto z=unit({b.x-a.x,b.y-a.y,b.z-a.z});
    const auto x=unit(cross(std::abs(z[2])<.85?V{0,0,1}:V{0,1,0},z));
    const auto y=cross(z,x);
    return {x[0],y[0],z[0],x[1],y[1],z[1],x[2],y[2],z[2]};
}
struct Channel {int m=-1;double fraction=0,centre=0;std::array<double,5> fractions{};};
Channel assess(const LocalAngularProjection& p,double determined){
    Channel out;out.centre=p.centre_mean_fraction;
    if(p.status!=MetricSubspaceStatus::Available || p.spins.size()!=1 ||
       !std::isfinite(p.angular_partition_residual) || std::abs(p.angular_partition_residual)>2e-5 ||
       !(p.centre_mean_fraction>=OrbitalChemistryOptions{}.pair_atom_weight_floor) ||
       p.represented_spin_orbital_rank!=1 || !(p.centre_projection_trace>1e-8))return out;
    for(int l=0;l<=4;++l)for(int c=0;c<2*l+1;++c){
        const auto value=p.component_projection_traces[l*l+c];
        if(std::isfinite(value))out.fractions[(c+1)/2]+=value/p.centre_projection_trace;
    }
    const auto best=std::max_element(out.fractions.begin(),out.fractions.end());
    out.fraction=*best;if(*best>=determined)out.m=static_cast<int>(best-out.fractions.begin());
    return out;
}
}
void annotate_nbo_bond_channels(NboIntegration& data,const Wavefunction& canonical){
    Wavefunction view=canonical;view.orbitals.clear();
    std::vector<std::size_t> descriptors;
    for(std::size_t i=0;i<data.orbitals.size();++i){const auto& o=data.orbitals[i];
        if(o.ref.kind!=NboOrbitalKind::NBO)continue;
        MolecularOrbital mo;mo.coefficients=o.coefficients;mo.spin=o.ref.spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;
        mo.occupation=o.occupation.value_or(0);view.orbitals.push_back(std::move(mo));descriptors.push_back(i);
    }
    using Key=std::tuple<std::size_t,std::size_t,std::size_t>;
    std::map<Key,std::unique_ptr<LocalAngularProjectionWorkspace>> workspaces;
    constexpr const char* names[]={"sigma","pi","delta","phi","gamma"};
    std::map<std::tuple<int,std::size_t>,NboChannelEvidence> channels;
    for(std::size_t i=0;i<descriptors.size();++i){auto& o=data.orbitals[descriptors[i]];
        if(o.atoms.size()!=2 || o.atoms[0]>=view.atoms.size() || o.atoms[1]>=view.atoms.size())continue;
        const auto pair=std::minmax(o.atoms[0],o.atoms[1]);
        const auto axes=frame(view.atoms[pair.first],view.atoms[pair.second]);
        std::array<Channel,2> parts;std::string error;
        for(std::size_t side=0;side<2;++side){const auto atom=side?pair.second:pair.first;const Key key{pair.first,pair.second,atom};
            try {auto& workspace=workspaces[key];if(!workspace)workspace=std::make_unique<LocalAngularProjectionWorkspace>(view,atom,axes);
                parts[side]=assess(workspace->project(std::span<const std::size_t>(&i,1)),OrbitalChemistryOptions{}.determined_fraction);
            }catch(const std::exception& e){error=e.what();}
        }
        const bool assigned=error.empty() && parts[0].m>=0 && parts[0].m==parts[1].m;
        const std::string label=assigned?names[parts[0].m]:"mixed/unresolved";
        NboChannelEvidence evidence;
        evidence.status=assigned?"available":"insufficient_evidence";
        evidence.channel=assigned?label:"unassigned";
        evidence.method="S-metric centre-conditioned bond-axis angular projection";
        evidence.reason=assigned?"Both centres support the same angular channel":
            error.empty()?"No common angular channel passes both centre thresholds":error;
        evidence.spin=o.ref.spin;
        evidence.axis_atoms={pair.first,pair.second};
        evidence.classification_threshold=OrbitalChemistryOptions{}.determined_fraction;
        evidence.centre_threshold=OrbitalChemistryOptions{}.pair_atom_weight_floor;
        evidence.source=o.source;
        for(std::size_t side=0;side<2;++side){
            evidence.centre_fractions[side]=parts[side].centre;
            evidence.conditional_angular_fractions[side]=parts[side].fractions;
        }
        o.channels.push_back(evidence);
        std::ostringstream detail;detail<<std::setprecision(6)<<"Derived local bond-axis channel: "<<label
            <<". Axis atoms "<<pair.first+1<<"-"<<pair.second+1<<". Separate S-metric centre projections:";
        for(std::size_t side=0;side<2;++side){detail<<" atom "<<(side?pair.second:pair.first)+1<<" centre="<<parts[side].centre<<"; conditional |m| fractions=";
            for(auto fraction:parts[side].fractions)detail<<fraction<<',';}
        detail<<". Derived label threshold="<<OrbitalChemistryOptions{}.determined_fraction
            <<"; required centre fraction="<<OrbitalChemistryOptions{}.pair_atom_weight_floor
            <<". Centre projections are not added as orthogonal populations. "<<error;
        o.detail+=(o.detail.empty()?"":"\n")+detail.str();
        if(assigned){o.display_label=o.label+" ["+label+"; derived]";
            channels[{static_cast<int>(o.ref.spin),o.ref.index}]=evidence;}
    }
    for(auto& evidence:data.structure){if(evidence.kind!="bond")continue;
        std::map<std::string,unsigned> counts;
        for(const auto& r:evidence.orbitals){if(r.kind!=NboOrbitalKind::NBO)continue;
            // Show per-orbital labels in the UI; aggregation never changes Lewis multiplicity.
            const auto found=channels.find({static_cast<int>(r.spin),r.index});if(found!=channels.end()){
                ++counts[found->second.channel];evidence.channels.push_back(found->second);
            }}
        if(!counts.empty()){std::string explanation="Derived channel labels on associated orbitals: ";
            for(const auto& [name,count]:counts)explanation+=name+" ("+std::to_string(count)+") ";
            evidence.detail+="\n"+explanation;}
    }
}
}
