#pragma once
#include "cov/nbo_aomo_labels.hpp"
#include <optional>

namespace cov::ui {
// An adapter for the existing typed field renderer, without claiming any NBO capability.
NboIntegration canonical_component_dataset(const Wavefunction&);
bool is_symmetry_component(const NboOrbitalSelection&);
std::optional<NboOrbitalSelection> symmetry_component_selection(
    const Wavefunction&,const NboIntegration*,const NboSalcModel*,
    const NboAomoName&,std::size_t component,
    std::optional<std::size_t> canonical_index={},
    std::optional<std::size_t> salc_index={});
} // namespace cov::ui
