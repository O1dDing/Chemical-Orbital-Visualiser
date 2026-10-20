#include "cov/mo_sigma_framework.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/local_angular_projection.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <tuple>

namespace cov {namespace {
using Vec=std::array<double,3>;
Vec cross(Vec a,Vec b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
Vec unit(Vec v){const double norm=std::hypot(v[0],v[1],v[2]);if(norm<1e-14)return {};for(auto& x:v)x/=norm;return v;}
std::array<double,9> bond_frame(const Atom& ligand,const Atom& centre) {
    const auto z=unit({centre.x-ligand.x,centre.y-ligand.y,centre.z-ligand.z});
    const auto x=unit(cross(std::abs(z[2])<.85?Vec{0,0,1}:Vec{0,1,0},z));
    const auto y=cross(z,x);return {x[0],y[0],z[0],x[1],y[1],z[1],x[2],y[2],z[2]};
}
const NboMatrix* transform(const NboDataset& dataset,NboSpin spin) {
    const NboMatrix* result=nullptr;
    for(const auto& m:dataset.matrices)if(m.kind=="NAONBO" && m.spin==spin) {
        if(result)return nullptr;result=&m;
    }
    return result;
}
}

MOSigmaFramework analyse_mo_sigma_framework(const Wavefunction& wf,
    const NboIntegration* integration,const RoutedAnalysis* route,
    const std::vector<std::size_t>& centres) {
    MOSigmaFramework out;
    if(!integration || !route || centres.empty() || !route->interaction_graph.available() ||
       integration->canonical_fingerprint!=route->canonical_fingerprint ||
       route->canonical_fingerprint!=nbo_canonical_fingerprint(wf) ||
       !integration->dataset.association.compatible) {
        out.detail="same-source-complete-NBO-or-centre-connectivity-unavailable";return out;
    }
    const auto n=static_cast<std::size_t>(wf.basis_count);
    if(!n || wf.ao_overlap.size()!=n*n) {
        out.detail="complete-AO-metric-unavailable";return out;
    }
    std::set<std::pair<std::size_t,std::size_t>> axes;
    for(const auto& edge:route->interaction_graph.value->edges)
        if(edge.strength==InteractionStrength::StrongConnectivity)for(const auto centre:centres) {
            const auto donor=edge.atom_a==centre?edge.atom_b:edge.atom_b==centre?edge.atom_a:wf.atoms.size();
            if(centre<wf.atoms.size() && donor<wf.atoms.size() &&
               std::find(centres.begin(),centres.end(),donor)==centres.end())axes.emplace(centre,donor);
        }
    if(axes.empty()){out.detail="first-shell-donor-axes-unavailable";return out;}
    const bool separate_spin=std::any_of(integration->dataset.orbitals.begin(),integration->dataset.orbitals.end(),
        [](const auto& o){return o.spin!=NboSpin::Total;});
    Wavefunction source=wf;source.orbitals.clear();
    std::vector<const NboOrbitalDescriptor*> candidates;
    std::vector<const NboOrbital*> literals;
    for(const auto& descriptor:integration->orbitals) {
        if(descriptor.ref.kind!=NboOrbitalKind::NBO || !descriptor.orthonormal_basis ||
           descriptor.coefficients.size()!=n || (separate_spin && descriptor.ref.spin==NboSpin::Total))continue;
        const auto literal=std::find_if(integration->dataset.orbitals.begin(),integration->dataset.orbitals.end(),
            [&](const auto& o){return o.id==descriptor.ref.index+1 && o.spin==descriptor.ref.spin;});
        if(literal==integration->dataset.orbitals.end() ||
           (literal->kind!="LP" && literal->kind!="BD" && literal->kind!="3C") ||
           literal->occupation<.5*(literal->spin==NboSpin::Total?2.:1.))continue;
        const bool eligible=std::any_of(axes.begin(),axes.end(),[&](const auto& axis) {
            const auto has=[&](auto atom){return std::find(descriptor.atoms.begin(),descriptor.atoms.end(),atom)!=descriptor.atoms.end();};
            return has(axis.second) && (literal->kind=="LP"?descriptor.atoms.size()==1:has(axis.first));
        });
        if(!eligible)continue;
        MolecularOrbital mo;mo.coefficients=descriptor.coefficients;
        mo.spin=descriptor.ref.spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;mo.occupation=static_cast<float>(literal->occupation);
        source.orbitals.push_back(std::move(mo));candidates.push_back(&descriptor);literals.push_back(&*literal);
    }
    using Axis=std::pair<std::size_t,std::size_t>;
    std::map<Axis,std::unique_ptr<LocalAngularProjectionWorkspace>> workspaces;
    for(std::size_t i=0;i<candidates.size();++i) {
        const auto& descriptor=*candidates[i];const auto& literal=*literals[i];
        const auto* matrix=transform(integration->dataset,descriptor.ref.spin);
        if(!matrix || matrix->rows!=n || matrix->columns<=descriptor.ref.index ||
           matrix->values.size()!=matrix->rows*matrix->columns)continue;
        double best=0;MOSigmaDonorSource selected;bool available=false;
        for(const auto& [centre,donor]:axes) {
            const auto has=[&](auto atom){return std::find(descriptor.atoms.begin(),descriptor.atoms.end(),atom)!=descriptor.atoms.end();};
            if(!has(donor) || (literal.kind!="LP" && !has(centre)))continue;
            double atom_weight=0,valence_sp=0;
            for(const auto& nao:integration->dataset.naos)if(nao.spin==descriptor.ref.spin && nao.atom==donor+1 && nao.id && nao.id<=n) {
                const double coefficient=matrix->values[(nao.id-1)*matrix->columns+descriptor.ref.index];
                const double weight=coefficient*coefficient;atom_weight+=weight;
                if(nao.type.rfind("Val",0)==0 && !nao.angular.empty() &&
                   (nao.angular.front()=='s' || nao.angular.front()=='p'))valence_sp+=weight;
            }
            // These are conditional source-type/angle classifications. They
            // never replace any full-canonical composition denominator.
            if(valence_sp<.05 || atom_weight<=0 || valence_sp/atom_weight<.5)continue;
            try {
                auto& workspace=workspaces[{centre,donor}];
                if(!workspace)workspace=std::make_unique<LocalAngularProjectionWorkspace>(source,donor,
                    bond_frame(wf.atoms[donor],wf.atoms[centre]));
                const auto projection=workspace->project(std::span<const std::size_t>(&i,1));
                if(projection.status!=MetricSubspaceStatus::Available ||
                   projection.represented_spin_orbital_rank!=1 || !std::isfinite(projection.angular_partition_residual) ||
                   std::abs(projection.angular_partition_residual)>2e-5)continue;
                const auto& p=projection.component_projection_traces;
                const double sp=p[0]+p[1]+p[2]+p[3];
                if(!std::isfinite(sp) || sp<=0)continue;
                const double sigma=(p[0]+p[1])/sp;
                if(sigma<.50 || sigma<=best)continue;
                best=sigma;available=true;selected={descriptor.ref,centre,donor,literal.kind,
                    literal.occupation,valence_sp,sigma};
            }catch(const std::exception&) {continue;}
        }
        if(available)out.sources.push_back(selected);
    }
    if(out.sources.empty()){out.detail="no-verified-first-shell-occupied-sigma-donors";return out;}
    // Complete selected NBO vectors form one fixed orthogonal projector. This
    // rejects duplicate/defective columns instead of normalizing them away.
    for(std::size_t i=0;i<out.sources.size();++i)for(std::size_t j=0;j<=i;++j) {
        const auto& a=out.sources[i].orbital;const auto& b=out.sources[j].orbital;
        if(a.spin!=b.spin)continue;const auto* m=transform(integration->dataset,a.spin);double dot=0;
        for(std::size_t row=0;row<n;++row)dot+=m->values[row*m->columns+a.index]*m->values[row*m->columns+b.index];
        out.source_orthogonality_error=std::max(out.source_orthogonality_error,std::abs(dot-(i==j?1.:0.)));
    }
    if(out.source_orthogonality_error>5e-5) {
        out.status="rejected";out.detail="source-sigma-NBO-projector-not-orthogonal";return out;
    }
    out.canonical_weights.assign(wf.orbitals.size(),0);
    // RO canonical coefficients are a verified shared spatial basis. The
    // integration's physical-spin link inventory need not include beta links
    // to its alpha-labelled spatial rows. Project the actual source columns
    // with the unchanged AO metric, using the declared occupation model.
    const bool shared=wf.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared;
    if(shared && separate_spin) {
        const auto has_spin=[&](NboSpin spin){return std::any_of(integration->dataset.orbitals.begin(),
            integration->dataset.orbitals.end(),[&](const auto& o){return o.spin==spin;});};
        if(!has_spin(NboSpin::Alpha) || !has_spin(NboSpin::Beta)) {
            out.detail="shared-spatial-source-requires-both-physical-spin-NBO-blocks";return out;
        }
        out.shared_spatial_average=true;
    }
    std::vector<std::vector<double>> metric_columns(wf.orbitals.size(),std::vector<double>(n));
    for(std::size_t index=0;index<wf.orbitals.size();++index) {
        if(wf.orbitals[index].coefficients.size()!=n) {
            out.detail="complete-canonical-AO-columns-unavailable";return out;
        }
        for(std::size_t row=0;row<n;++row)for(std::size_t col=0;col<n;++col)
            metric_columns[index][row]+=wf.ao_overlap[row*n+col]*wf.orbitals[index].coefficients[col];
    }
    std::map<std::tuple<int,std::size_t,std::size_t>,double> unique_links;
    std::set<NboSpin> source_spins;
    for(const auto& selected:out.sources) {
        source_spins.insert(selected.orbital.spin);
        const auto descriptor=std::find_if(integration->orbitals.begin(),integration->orbitals.end(),
            [&](const auto& d){return d.ref==selected.orbital;});
        for(std::size_t index=0;index<wf.orbitals.size();++index) {
            if(!shared && (selected.orbital.spin==NboSpin::Beta)!=(wf.orbitals[index].spin==Spin::Beta))continue;
            double coefficient=0;
            for(std::size_t row=0;row<n;++row)coefficient+=descriptor->coefficients[row]*metric_columns[index][row];
            unique_links[{static_cast<int>(selected.orbital.spin),selected.orbital.index,index}]=coefficient*coefficient;
        }
    }
    if(unique_links.empty()){out.detail="complete-NBO-canonical-transform-unavailable";return out;}
    for(const auto& selected:out.sources) {
        double norm=0;
        for(const auto& [key,weight]:unique_links)
            if(std::get<0>(key)==static_cast<int>(selected.orbital.spin) &&
               std::get<1>(key)==selected.orbital.index)norm+=weight;
        out.source_canonical_closure_error=std::max(out.source_canonical_closure_error,std::abs(norm-1));
        if(!std::isfinite(norm) || std::abs(norm-1)>5e-5) {
            out.detail="selected-sigma-NBO-column-not-complete-in-canonical-source-space";return out;
        }
    }
    const double spin_average=out.shared_spatial_average?.5:1.;
    for(const auto& [key,weight]:unique_links)out.canonical_weights[std::get<2>(key)]+=spin_average*weight;
    for(std::size_t i=0;i<out.canonical_weights.size();++i) {
        if(out.canonical_weights[i]>1+5e-5) {
            out.status="rejected";out.detail="sigma-projector-exceeds-original-canonical-norm";return out;
        }
        if(wf.orbitals[i].occupation>.05)out.occupied_trace+=out.canonical_weights[i];
    }
    for(const auto& source:out.sources) {
        if(source.orbital.spin==NboSpin::Alpha)++out.alpha_source_rank;
        else if(source.orbital.spin==NboSpin::Beta)++out.beta_source_rank;
        else ++out.total_source_rank;
    }
    out.rank=out.shared_spatial_average?std::max(out.alpha_source_rank,out.beta_source_rank):out.sources.size();
    out.available=true;out.status="available";
    out.detail="complete-orthogonal-occupied-LP-and-centre-ligand-BD-3C; actual-axis-Val-sp-sigma-classification";
    if(out.shared_spatial_average)out.detail+="; shared-spatial-alpha-beta-projector-average; reading-rank=max-physical-spin-source-ranks";
    return out;
}
}
