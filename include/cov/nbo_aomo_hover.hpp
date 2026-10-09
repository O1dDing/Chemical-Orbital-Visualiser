#pragma once

#include "cov/ui.hpp"
#include <string>
#include <vector>

namespace cov {
struct Wavefunction;
struct NboIntegration;
struct NboSalcModel;
namespace ui {
struct NboAomoNode;
// Presentation of existing evidence only. No chemistry is derived, and no
// orbital/selection is modified. Common spatial nodes also expose both source
// spin expectations; complete rotated source-member values are in the detail.
std::vector<std::string> nbo_aomo_hover_lines(const NboAomoNode&,
    const Wavefunction&, const NboIntegration&, const NboSalcModel*, Language);
std::string nbo_aomo_hover_glyph_seed(Language);
} // namespace ui
} // namespace cov
