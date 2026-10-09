#pragma once

#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <set>
#include <string>

namespace cov {

// Chemical applicability and numerical completeness answer separate questions.
// Full NAO data must never turn a non-metal centre into a metal centre.
struct CompositionDisplayPolicy {
    bool complete_nao = false;
    bool complete_metal_composition = false;
    bool metal_ligand = false;
    bool reference_populations = false;
    bool reference_shell_traces = false;
};

[[nodiscard]] inline CompositionDisplayPolicy composition_display_policy(
    const MOGroupCompositionLedger& composition,
    const MetalLigandDetailAvailability& scope) noexcept {
    CompositionDisplayPolicy result;
    result.complete_nao = composition.available && composition.complete;
    result.complete_metal_composition = result.complete_nao && scope.composition;
    result.metal_ligand = result.complete_metal_composition ||
        scope.scope == ChemistryStatus::Determined;
    result.reference_populations = scope.scope == ChemistryStatus::Determined &&
        !result.complete_nao && scope.populations;
    result.reference_shell_traces = !result.complete_nao;
    return result;
}

// The weight remains the complete-MO trace for this angular family. Multiple
// centres may contribute, so name every actual radial shell rather than
// assigning one element's principal quantum number to the entire sum.
[[nodiscard]] inline std::string current_composition_shell_caption(
    const Wavefunction& wavefunction,
    const std::vector<MOCurrentRadialShell>& shells, const int angular_l) {
    constexpr const char* families[] = {"s", "p", "d", "f"};
    if(angular_l<0 || angular_l>3)return {};
    std::string result;
    std::set<std::pair<std::size_t,int>> included;
    for(const auto& shell:shells) {
        if(shell.l!=angular_l || shell.n<=0 || shell.atom>=wavefunction.atoms.size() ||
           !included.emplace(shell.atom,shell.n).second)continue;
        if(!result.empty())result+=" + ";
        const auto& atom=wavefunction.atoms[shell.atom];
        result+=(atom.symbol.empty()?"Atom "+std::to_string(shell.atom+1):atom.symbol);
        result+=' '+std::to_string(shell.n)+families[angular_l];
    }
    return result.empty()?families[angular_l]:result;
}

} // namespace cov
