#pragma once
#include "cov/nbo_integration.hpp"
namespace cov {
// Derived, conditional bond-axis labels from actual S-metric projections.
// Producer labels and all orbital coefficients remain unchanged.
void annotate_nbo_bond_channels(NboIntegration&,const Wavefunction&);
}
