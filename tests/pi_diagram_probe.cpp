#include "cov/chemistry_route.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/nbo_salc.hpp"
#include "cov/wavefunction_io.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

// Small observation adapter for saved inputs; acceptance compares the actual
// channel/member records, not a process exit or a molecule-name lookup.
int main(int argc,char** argv) {try {
    if(argc!=4&&argc!=5)throw std::runtime_error("Usage: cov_pi_diagram_probe canonical.fchk analysis_directory|- output.json [without-e2|duplicate-e2|scope-only]");
    const auto w=cov::parse_wavefunction(argv[1]);
    const auto identity=cov::nbo_canonical_fingerprint(w);
    std::optional<cov::NboIntegration> integration;
    if(std::string(argv[2])!="-") {
        const auto found=cov::discover_nbo_inputs({std::filesystem::path(argv[2])});
        if(found.candidates.size()!=1)throw std::runtime_error("Exactly one NBO analysis is required");
        integration=cov::read_nbo_integration(w,found.candidates.front());
        cov::annotate_nbo_bond_channels(*integration,w);
        if(argc==5) {
            const std::string mode=argv[4];
            if(mode=="without-e2")integration->dataset.e2.clear();
            else if(mode=="duplicate-e2") {
                const auto original=integration->dataset.e2;
                integration->dataset.e2.insert(integration->dataset.e2.end(),original.begin(),original.end());
            }else if(mode!="scope-only")throw std::runtime_error("Unknown evidence fault-injection mode");
        }
    }
    std::optional<cov::NboSalcModel> model;
    if(integration)model=cov::build_nbo_salc_model(w,*integration);
    const auto route=cov::route_chemistry(w,integration?&*integration:nullptr);
    std::ofstream out(argv[3]);
    if(!out)throw std::runtime_error("Cannot write output");
    out<<"{\"route\":"<<cov::serialize_routed_analysis_json(route)<<",\"views\":[";
    for(int compact=0;compact<2;++compact) {
        cov::MODiagramOptions options;options.routed=&route;
        options.nbo_source=integration?&*integration:nullptr;
        options.routed_identity=route.canonical_fingerprint;
        options.hide_ligand_centred_intermediates=compact!=0;
        options.aomo_scope=compact?1:2;
        if(model)for(const auto& fragment:model->fragments)
            if(fragment.side==0&&fragment.atoms.size()==1)options.display_centre_atoms.push_back(fragment.atoms.front());
        if(options.display_centre_atoms.size()!=1)options.display_centre_atoms.clear();
        const auto diagram=cov::build_mo_diagram_data(w,options);
        if(compact)out<<',';
        out<<"{\"compact\":"<<(compact?"true":"false")
           <<",\"level_group_count\":"<<diagram.levels.size()
           <<",\"pi_interactions\":"<<((argc==5&&std::string(argv[4])=="scope-only")?"[]":cov::orbital_energy_gap_array_json(diagram.pi_interactions))
           <<",\"candidates\":"<<((argc==5&&std::string(argv[4])=="scope-only")?"[]":cov::pi_partner_candidates_json(diagram.pi_partner_candidates))
           <<",\"pi_mode_networks\":[";
        for(std::size_t i=0;i<diagram.pi_mode_networks.size();++i){if(i)out<<',';
            out<<cov::pi_mode_network_assessment_json(diagram.pi_mode_networks[i]);}
        out<<"],\"visible_members\":[";
        bool comma=false;
        for(const auto& row:diagram.levels)for(const auto member:row.member_indices) {
            if(comma)out<<',';comma=true;out<<member;
        }
        out<<"],\"scoped_frozen\":[";bool first_channel=true;
        for(const auto& result:route.pi_couplings)if(result.available()){
            const auto& channel=*result.value;std::vector<std::size_t> domain;
            for(const auto& level:diagram.levels)for(auto member:level.member_indices)
                if(channel.canonical_members_are_verified_shared_spatial?w.orbitals[member].spin==cov::Spin::Alpha:
                   ((channel.spin==cov::NboSpin::Beta)==(w.orbitals[member].spin==cov::Spin::Beta)))domain.push_back(member);
            cov::PiPartnerAssessment evidence;evidence.channel.channel_id=channel.id;
            evidence.channel.frozen_operator=cov::assess_pi_frozen_display_scope(channel,domain);
            if(!first_channel)out<<',';first_channel=false;out<<cov::pi_partner_assessment_json(evidence);
        }
        out<<"]}";
    }
    out<<"],\"canonical_immutable\":"<<(identity==cov::nbo_canonical_fingerprint(w)?"true":"false")<<'}';
    return out?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
