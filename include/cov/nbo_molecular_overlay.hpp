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
#include "cov/molecular_overlay.hpp"
#include "cov/nbo_integration.hpp"
#include "cov/interaction_graph.hpp"
#include "cov/chemistry_route.hpp"
namespace cov {
struct NboElectronicSymmetryScope;
MoleculeOverlay make_nbo_molecule_overlay(const NboIntegration&,const InteractionGraph&,
    std::size_t atom_count,const std::vector<std::size_t>& selected_atoms,
    std::optional<std::size_t> selected_structure,AtomScalarMode,bool show_indices,bool show_e2,
    const RoutedAnalysis* routed=nullptr,const Wavefunction* canonical=nullptr,
    NboBondDisplayMode bond_mode=NboBondDisplayMode::DefaultSkeleton,
    const NboElectronicSymmetryScope* electronic_scope=nullptr);
std::string serialize_molecule_overlay_scalars_json(const MoleculeOverlay&,
                                                     const RoutedAnalysis* routed=nullptr);
}
