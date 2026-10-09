#include "cov/chemistry_route.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/nbo_molecular_overlay.hpp"
#include "cov/wavefunction_io.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
std::string json_string(const std::string& value){
    std::ostringstream out;out<<'"';
    for(unsigned char c:value){
        if(c=='"'||c=='\\')out<<'\\'<<c;
        else if(c=='\n')out<<"\\n";
        else if(c=='\r')out<<"\\r";
        else if(c=='\t')out<<"\\t";
        else if(c<0x20){
            constexpr char hex[]="0123456789abcdef";
            out<<"\\u00"<<hex[c>>4]<<hex[c&15];
        }else out<<c;
    }
    return out.str()+'"';
}
}

// Observation adapter: this calls exactly the producer-independent router
// used by the application. An independent raw-input oracle judges its values.
int main(int argc,char** argv) {try {
    if(argc!=4)throw std::runtime_error(
        "Usage: cov_chemistry_route_probe canonical.fchk analysis_directory|- output.json");
    const auto canonical=cov::parse_wavefunction(argv[1]);
    const auto before=cov::nbo_canonical_fingerprint(canonical);
    std::optional<cov::NboIntegration> integration;
    if(std::string(argv[2])!="-") {
        const auto found=cov::discover_nbo_inputs({std::filesystem::path(argv[2])});
        if(found.candidates.size()!=1)
            throw std::runtime_error("Exactly one analysis candidate is required");
        integration=cov::read_nbo_integration(canonical,found.candidates.front());
        cov::annotate_nbo_bond_channels(*integration,canonical);
    }
    const auto route=cov::route_chemistry(canonical,integration?&*integration:nullptr);
    std::string nbo_render_status=integration?"unavailable":"not_analysed";
    std::string nbo_render_error;
    std::size_t nbo_alpha=0,nbo_beta=0;
    if(integration){
        const auto* capability=cov::nbo_capability(*integration,"nbo");
        if(capability&&capability->available()){
            try{
                const auto rendered=cov::make_nbo_wavefunction(integration->dataset,canonical);
                for(const auto& orbital:rendered.orbitals)
                    if(orbital.spin==cov::Spin::Beta)++nbo_beta;else ++nbo_alpha;
                nbo_render_status="available";
            }catch(const std::exception& error){
                nbo_render_status="rejected";nbo_render_error=error.what();
            }
        }
        else if(capability){
            nbo_render_status=cov::nbo_capability_state_name(capability->state);
            nbo_render_error=capability->detail;
        }
    }
    cov::ActiveOrbitalView active;
    active.kind=cov::ActiveOrbitalKind::Canonical;
    active.source_id=route.canonical_fingerprint;
    active.semantic_kind="canonical";
    if(!canonical.orbitals.empty()) {
        active.canonical_index=0;active.rendered_index=0;
        active.label="MO 1";
        active.source_label=active.label;
        active.display_name_evidence="canonical source ordinal; no presentation override";
        active.spin=canonical.orbital_occupation_model==cov::OrbitalOccupationModel::CanonicalShared?
            cov::NboSpin::Total:canonical.orbitals.front().spin==cov::Spin::Beta?
            cov::NboSpin::Beta:cov::NboSpin::Alpha;
        active.source_spin=canonical.orbitals.front().spin==cov::Spin::Beta?
            cov::NboSpin::Beta:cov::NboSpin::Alpha;
        active.spin_semantics=active.spin==cov::NboSpin::Total?
            "shared spatial-orbital display; alpha is the producer channel field":
            "explicit canonical spin orbital";
    }
    const auto output=std::filesystem::path(argv[3]);
    if(output.has_parent_path())std::filesystem::create_directories(output.parent_path());
    std::ofstream out(output,std::ios::binary);
    if(!out)throw std::runtime_error("Cannot write routed observation");
    out<<"{\"schema\":\"cov.chemistry.route-probe.v1\",\"route\":"
       <<cov::serialize_routed_analysis_json(route,true)
       <<",\"active_view\":"<<cov::serialize_active_orbital_view_json(active)
       <<",\"integration\":";
    if(integration)out<<cov::serialize_nbo_integration_json(*integration);
    else out<<"null";
    out<<",\"nbo_render\":{\"status\":"<<json_string(nbo_render_status)
       <<",\"error\":"<<(nbo_render_error.empty()?"null":json_string(nbo_render_error))
       <<",\"alpha_orbitals\":"<<nbo_alpha
       <<",\"beta_orbitals\":"<<nbo_beta
       <<",\"total_orbitals\":"<<(nbo_alpha+nbo_beta)
       <<",\"recheck_scope\":\"matched archive plus the AONBO matrix family\"}";
    if(integration) {
        const auto overlay=cov::make_nbo_molecule_overlay(*integration,*route.interaction_graph.value,
            canonical.atoms.size(),{},std::nullopt,cov::AtomScalarMode::Element,false,false,&route);
        out<<",\"molecule_overlay\":{\"bonds\":[";
        bool comma=false;
        for(const auto& bond:overlay.bonds){
            if(comma)out<<',';comma=true;
            out<<"{\"atoms\":["<<bond.atom_a<<','<<bond.atom_b<<"],\"style\":"
               <<static_cast<int>(bond.style)<<",\"multiplicity\":"<<bond.multiplicity<<'}';
        }
        out<<"],\"multicentre\":[";comma=false;
        for(const auto& group:overlay.multicentre){
            if(comma)out<<',';comma=true;out<<"{\"atoms\":[";
            for(std::size_t i=0;i<group.atoms.size();++i){if(i)out<<',';out<<group.atoms[i];}
            out<<"],\"supporting_records\":"<<group.evidence_indices.size()<<'}';
        }
        out<<"]}";
    }
    out<<",\"canonical_immutable\":"
       <<(before==cov::nbo_canonical_fingerprint(canonical)?"true":"false")<<'}';
    if(!out)throw std::runtime_error("Routed observation write failed");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
