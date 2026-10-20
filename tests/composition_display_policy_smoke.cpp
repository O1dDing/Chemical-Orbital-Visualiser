#include "cov/composition_display_policy.hpp"

#include <iostream>
#include <stdexcept>

namespace {
void require(const bool condition,const char* message) {
    if(!condition)throw std::runtime_error(message);
}
}

int main() {try {
    cov::MOGroupCompositionLedger complete;
    complete.available=complete.complete=true;
    cov::MetalLigandDetailAvailability applicability;
    applicability.populations=true;

    // A complete numerical source cannot authorize chemical terminology:
    // benzene, water and PPh3 have complete NAO data but no metal scope.
    for(const auto scope:{cov::ChemistryStatus::NotApplicable,cov::ChemistryStatus::Unavailable}) {
        applicability.scope=scope;
        const auto p=cov::composition_display_policy(complete,applicability);
        require(p.complete_nao && !p.complete_metal_composition && !p.metal_ligand && !p.reference_populations,
            "Complete NAO source bypasses chemical applicability");
        require(!p.reference_shell_traces,
            "Complete NAO source retains conflicting local-reference percentages");
    }

    // A generic metal centre (including alkali metals) with an associated
    // ligand scope does not require a d-block local point-group analysis.
    applicability.composition=true;
    applicability.scope=cov::ChemistryStatus::NotApplicable;
    auto p=cov::composition_display_policy(complete,applicability);
    require(p.metal_ligand && p.complete_metal_composition && p.complete_nao && !p.reference_populations,
        "Metal composition incorrectly requires d-block ligand-field applicability");
    auto generic_partial=complete;generic_partial.complete=false;
    p=cov::composition_display_policy(generic_partial,applicability);
    require(!p.metal_ligand && !p.reference_populations,
        "Partial generic metal composition borrows unavailable local reference data");
    applicability.composition=false;

    // Chemical scope may remain valid when detailed interaction quantities
    // are absent. Quantitative source choice is independent of that scope.
    applicability.scope=cov::ChemistryStatus::Determined;
    applicability.composition=true;
    applicability.populations=false;
    p=cov::composition_display_policy(complete,applicability);
    require(p.metal_ligand && p.complete_metal_composition && p.complete_nao && !p.reference_populations,
        "Complete metal composition requires unrelated reference populations");
    applicability.populations=true;
    for(const auto fractions:{.501492,.497527,.400894,.556189}) {
        complete.centre_current_d=fractions;
        const auto result=cov::composition_display_policy(complete,applicability);
        require(result.metal_ligand && result.complete_metal_composition && result.complete_nao && !result.reference_populations &&
                !result.reference_shell_traces,
            "Composition fraction changes source or chemical display policy");
    }
    cov::MOGroupCompositionLedger partial=complete;
    partial.complete=false;
    p=cov::composition_display_policy(partial,applicability);
    require(p.metal_ligand && !p.complete_nao && p.reference_populations && p.reference_shell_traces,
        "Partial data silently presented as complete NAO composition");
    applicability.populations=false;
    p=cov::composition_display_policy(partial,applicability);
    require(!p.reference_populations,"Unavailable reference becomes a measured zero");

    // A current 4s framework and current 3d role must remain distinguishable,
    // independent of view numbering. Multi-centre sums retain every identity.
    cov::Wavefunction w;
    w.atoms.resize(2);
    w.atoms[0].symbol="Fe";w.atoms[0].atomic_number=26;
    w.atoms[1].symbol="Co";w.atoms[1].atomic_number=27;
    const std::vector<cov::MOCurrentRadialShell> shells={
        {0,4,0,{}},{0,4,1,{}},{0,3,2,{}},{0,3,2,{}},
        {1,3,2,{}},{99,9,2,{}},{0,0,2,{}}};
    require(cov::current_composition_shell_caption(w,shells,0)=="Fe 4s",
        "Current s contribution loses principal radial-shell identity");
    require(cov::current_composition_shell_caption(w,shells,1)=="Fe 4p",
        "Current p contribution loses principal radial-shell identity");
    require(cov::current_composition_shell_caption(w,shells,2)=="Fe 3d + Co 3d",
        "Multi-centre d sum is assigned to only one centre or includes malformed shells");
    require(cov::current_composition_shell_caption(w,{},2)=="d",
        "Missing radial identity invents a principal quantum number");
    std::cout<<"Composition applicability, exclusive numerical source, partial-data and radial-shell display checks passed\n";
    return 0;
} catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';return 1;
}}
