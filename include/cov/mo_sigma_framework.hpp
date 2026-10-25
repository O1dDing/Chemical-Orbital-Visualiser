/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

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
