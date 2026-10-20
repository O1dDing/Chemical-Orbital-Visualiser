#pragma once
#include "cov/nbo_integration.hpp"
#include <string>
#include <vector>
namespace cov {
struct RoutedAnalysis;
struct MOSigmaDonorSource {
    NboOrbitalRef orbital;
    std::size_t centre=0, ligand=0;
    std::string kind;
    double occupation=0, ligand_valence_sp_weight=0, axis_sigma_fraction=0;
};
// An overlapping channel projector, not an extra composition bucket. Its
// weights always retain the full original canonical norm and complete NBOs.
struct MOSigmaFramework {
    bool available=false;
    std::string status="unavailable", detail;
    std::string source="complete-same-source-occupied-NBO-sigma-projector";
    std::size_t rank=0;
    std::size_t alpha_source_rank=0, beta_source_rank=0, total_source_rank=0;
    bool shared_spatial_average=false;
    double source_orthogonality_error=0, source_canonical_closure_error=0, occupied_trace=0;
    double retained_occupied_trace=0, trace_target=.90;
    std::size_t display_member_budget=0, retained_members=0;
    std::vector<MOSigmaDonorSource> sources;
    std::vector<double> canonical_weights;
};
MOSigmaFramework analyse_mo_sigma_framework(const Wavefunction&,
    const NboIntegration*,const RoutedAnalysis*,const std::vector<std::size_t>& centres);
}
