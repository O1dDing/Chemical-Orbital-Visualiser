#include "cov/chemistry_route.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/mo_group_display_json.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/nbo_salc.hpp"
#include "cov/wavefunction_io.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <stdexcept>

namespace {
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
double current(const cov::MOGroupCompositionLedger& c) {
    return c.centre_current_s+c.centre_current_p+c.centre_current_d+c.centre_current_f;
}
bool paired(const cov::MODiagramGroupAudit& g) {
    return std::find(g.display_decision.reason_codes.begin(),g.display_decision.reason_codes.end(),
        "paired-opposite-spin-source-group")!=g.display_decision.reason_codes.end();
}
void indices(std::ostream& out,const std::vector<std::size_t>& values) {
    out<<'[';for(std::size_t i=0;i<values.size();++i){if(i)out<<',';out<<values[i];}out<<']';
}
}

// Observe saved wavefunctions through the actual production builder. Curves
// describe source norm/trace and reading budgets; they are not energy estimates.
int main(int argc,char** argv) {try {
    require(argc==4,"Usage: cov_mo_brief_scope_probe canonical.fchk analysis_directory output.json");
    const auto wf=cov::parse_wavefunction(argv[1]);
    const auto identity=cov::nbo_canonical_fingerprint(wf);
    const auto found=cov::discover_nbo_inputs({argv[2]});
    require(found.candidates.size()==1,"Exactly one associated NBO analysis required");
    auto integration=cov::read_nbo_integration(wf,found.candidates.front());
    cov::annotate_nbo_bond_channels(integration,wf);
    const auto model=cov::build_nbo_salc_model(wf,integration);
    const auto route=cov::route_chemistry(wf,&integration);
    cov::MODiagramOptions options;options.routed=&route;options.nbo_source=&integration;
    options.hide_ligand_centred_intermediates=true;
    for(const auto& fragment:model.fragments)
        if(fragment.side==0 && fragment.atoms.size()==1)options.display_centre_atoms.push_back(fragment.atoms.front());
    options.selected_index=wf.orbitals.empty()?0:wf.orbitals.size()-1;
    std::ofstream out(argv[3]);require(bool(out),"Cannot open output");out<<std::setprecision(17);
    out<<"{\"canonical_fingerprint\":"<<cov::mo_display_json_detail::quote(identity)
       <<",\"policy\":{\"version\":\"current-shell-reading-budget-v1\",\"minimum_current_coverage\":0.05,"
         "\"energy_kernel_hartree\":0.75,\"weighted_trace_target\":0.90,"
         "\"virtual_members_per_angular_family_per_centre\":\"2*(2*l+1)\","
         "\"semantics\":\"finite low-energy reading budget; complete groups; major relations separately protected\"},\"modes\":[";
    std::set<std::size_t> brief_members;
    cov::MODiagramData brief;
    for(unsigned scope=1;scope<=3;++scope) {
        options.aomo_scope=scope;
        const auto data=cov::build_mo_diagram_data(wf,options);
        require(data.selection.counts_are_final,"Counts are not from final graph");
        std::set<std::size_t> members;
        for(const auto& level:data.levels) {
            members.insert(level.member_indices.begin(),level.member_indices.end());
            for(auto counterpart:level.member_spin_counterparts)
                if(counterpart<wf.orbitals.size())members.insert(counterpart);
        }
        require(members.size()==data.selection.final_member_count,"Final member count differs from graph");
        require(data.levels.size()==data.selection.final_group_count,"Final group count differs from graph");
        require(data.selection.valence_occupied_count+data.selection.frontier_virtual_count==members.size(),
            "Occupied/empty final counts do not close");
        if(scope==3)require(members.size()==wf.orbitals.size(),"Full mode did not restore all canonical source members");
        if(scope==1){brief=data;brief_members=members;}
        if(scope>1)out<<',';
        out<<"{\"scope\":"<<scope<<",\"group_count\":"<<data.levels.size()
           <<",\"member_count\":"<<members.size()<<",\"occupied_members\":"<<data.selection.valence_occupied_count
           <<",\"empty_members\":"<<data.selection.frontier_virtual_count<<",\"members\":";
        indices(out,data.selection.included_indices);
        out<<",\"groups\":[";
        for(std::size_t i=0;i<data.group_audit.size();++i) {
            if(i)out<<',';const auto& group=data.group_audit[i];
            if(group.composition.complete) {
                const auto& c=group.composition;
                require(std::abs(current(c)+c.centre_other+c.ligand_valence+c.ligand_other+c.other_atoms+c.core+c.unresolved-1)<1e-5,
                    "Complete source ledger does not close");
            }
            out<<"{\"members\":";indices(out,group.member_indices);
            out<<",\"energy_hartree\":"<<group.energy_hartree<<",\"total_occupation\":"<<group.total_occupation
               <<",\"composition\":"<<cov::mo_group_composition_json(group.composition)
               <<",\"display_decision\":"<<cov::mo_group_display_decision_json(group.display_decision)<<'}';
        }
        out<<"],\"sigma_framework\":"<<cov::mo_sigma_framework_json(data.sigma_framework)
           <<",\"pi_mode_networks\":[";
        for(std::size_t i=0;i<data.pi_mode_networks.size();++i) {
            if(i)out<<',';out<<cov::pi_mode_network_assessment_json(data.pi_mode_networks[i]);
        }
        out<<"],\"pair_count\":"<<data.pi_interactions.size()<<'}';
    }
    options.aomo_scope=1;options.selected_index=0;
    const auto reselected=cov::build_mo_diagram_data(wf,options);
    require(std::set<std::size_t>(reselected.selection.included_indices.begin(),reselected.selection.included_indices.end())==brief_members,
        "Brief source membership depends on selection");
    out<<"],\"coverage_curve\":[";
    bool comma=false;
    for(double threshold:{0.0,0.01,0.025,0.05,0.10,0.15,0.20,0.30}) {
        double total=0,retained=0;std::size_t groups=0,members=0;
        for(const auto& group:brief.group_audit)if(!paired(group) && group.composition.complete) {
            const auto& c=group.composition;const double trace=current(c)*c.member_count;
            total+=trace;
            if(current(c)>=threshold){retained+=trace;++groups;members+=c.member_count;}
        }
        if(comma)out<<',';comma=true;
        out<<"{\"threshold\":"<<threshold<<",\"groups\":"<<groups<<",\"members\":"<<members
           <<",\"current_trace\":"<<retained<<",\"total_current_trace\":"<<total
           <<",\"lost_fraction\":"<<(total>0?1-retained/total:0)<<'}';
    }
    double effective=std::numeric_limits<double>::infinity();
    for(const auto& group:brief.group_audit)if(!paired(group) && group.composition.complete &&
        group.total_occupation<=options.filter.occupation_threshold && current(group.composition)>=.05)
        effective=std::min(effective,group.energy_hartree);
    out<<"],\"low_energy_recovery_curve\":[";
    for(int l=0;l<=3;++l) {
        if(l)out<<',';double total=0,retained=0,raw=0,kept_raw=0;
        for(const auto& group:brief.group_audit)if(!paired(group) && group.composition.complete &&
            group.total_occupation<=options.filter.occupation_threshold && std::isfinite(effective)) {
            const auto& c=group.composition;
            const double w=l==0?c.centre_current_s:l==1?c.centre_current_p:l==2?c.centre_current_d:c.centre_current_f;
            const double trace=w*c.member_count;
            const double gap=std::max(0.0,group.energy_hartree-effective);
            const double weighted=trace/(1+std::pow(gap/.75,2));
            total+=weighted;raw+=trace;
            if(group.display_decision.included){retained+=weighted;kept_raw+=trace;}
        }
        out<<"{\"l\":"<<l<<",\"full_virtual_trace\":"<<raw<<",\"retained_virtual_trace\":"<<kept_raw
           <<",\"low_energy_weighted_trace\":"<<total<<",\"retained_weighted_trace\":"<<retained
           <<",\"weighted_recovery\":"<<(total>0?retained/total:1)<<'}';
    }
    out<<"],\"frozen_display_scopes\":[";comma=false;
    for(const auto& entry:route.pi_couplings)if(entry.available()) {
        const auto& channel=*entry.value;std::vector<std::size_t> domain;
        for(auto member:brief.selection.included_indices)
            if(channel.canonical_members_are_verified_shared_spatial?wf.orbitals[member].spin==cov::Spin::Alpha:
               ((channel.spin==cov::NboSpin::Beta)==(wf.orbitals[member].spin==cov::Spin::Beta)))domain.push_back(member);
        cov::PiPartnerAssessment evidence;evidence.channel.channel_id=channel.id;
        evidence.channel.frozen_operator=cov::assess_pi_frozen_display_scope(channel,domain);
        if(comma)out<<',';comma=true;out<<cov::pi_partner_assessment_json(evidence);
    }
    require(identity==cov::nbo_canonical_fingerprint(wf),"Canonical source changed during display");
    out<<"],\"canonical_immutable\":true,\"selection_independent\":true}";
    return out?0:2;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
