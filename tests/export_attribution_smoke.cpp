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

// This fixture deliberately uses only pre-existing APIs so it can also be
// linked against the release libraries for a byte/pixel regression comparison.
#include "cov/mo_diagram.hpp"
#include "cov/nbo_ui.hpp"
#include "cov/nbo_aomo_ui.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        const auto root = argc > 1 ? std::filesystem::path(argv[1]) :
            std::filesystem::current_path() / "export-attribution-artifacts";
        std::filesystem::create_directories(root);
        cov::MODiagramData mo;
        mo.selection.summary = "Synthetic attribution regression fixture";
        mo.selection.included_indices = {0, 1};
        for (std::size_t i = 0; i < 2; ++i) {
            cov::MODiagramLevel level;
            level.metadata.orbital_index = i;
            level.metadata.raw_mo_number = i + 1;
            level.metadata.display_label = std::to_string(i + 1);
            level.metadata.energy_hartree = i ? .10 : -.20;
            level.metadata.occupation = i ? 0.f : 2.f;
            level.metadata.symmetry = "A1";
            level.layout_energy_hartree = level.metadata.energy_hartree;
            level.member_indices = {i};
            level.electrons = i ? cov::ElectronGlyphs{} : cov::ElectronGlyphs{1, 1};
            mo.metadata.push_back(level.metadata);
            mo.annotations.push_back(level.annotation);
            mo.levels.push_back(level);
        }
        cov::MODiagramOptions options;
        options.width = 760; options.height = 620;
        mo.energy_transform = cov::build_energy_transform(
            {-.20, .10}, cov::EnergyAxisMode::NonlinearFocus, .055);
        std::string error;
        if (!cov::write_mo_diagram_svg(mo, options, root / "fixture.mo.svg", &error) ||
            !cov::write_mo_diagram_png(mo, options, root / "fixture.mo.png", &error))
            throw std::runtime_error(error);

        const auto diagram = cov::make_mo_diagram_view_snapshot(mo, options, 0);
        cov::NboDataset nbo;
        nbo.source.path = "PRIVATE_INPUT_PATH_SENTINEL";
        cov::NboOrbital orbital;
        orbital.id = 1; orbital.ordinal = 1; orbital.kind = "BD";
        orbital.label = "BD (1) H 1 - H 2"; orbital.occupation = 1.98;
        nbo.orbitals.push_back(orbital);
        cov::NboMoDecomposition decomposition;
        decomposition.available = true; decomposition.status = "available";
        cov::NboNaoContribution row;
        row.nao_id = 1; row.atom = 1; row.symbol = "H"; row.type = "Val";
        row.angular = "1s"; row.principal_n = 1; row.angular_l = 0;
        row.coefficient = 1.; row.weight = 1.;
        decomposition.rows.push_back(row);
        nbo.mo_decompositions.push_back(decomposition);
        cov::ui::NboFocusUIState focus;
        focus.canonical_index = 0; focus.visible_atoms.insert(1);
        cov::ActiveOrbitalView active;
        active.label = "MO 1"; active.canonical_index = 0; active.rendered_index = 0;
        cov::ui::export_nbo_bundle(nbo, 0, root / "fixture", nullptr, 0, .01,
            active, nullptr, &diagram, &focus);

        cov::NboIntegration integration;
        integration.id = "synthetic-attribution-fixture";
        cov::ui::NboAomoViewSnapshot view;
        view.integration_id = integration.id; view.capability_status = "available";
        view.canvas_width = 760; view.canvas_height = 480; view.footer_y = 450;
        view.paper_export = true;
        view.captions.push_back({"title", "AO / MO attribution fixture", 24, 15, 700, 30});
        for (int i = 0; i < 2; ++i) {
            cov::ui::NboAomoNode node;
            node.id = "fixture-node-" + std::to_string(i);
            node.label = node.individual_label = i ? "MO 1" : "H1 1s";
            node.label_atom_indices = i == 0;
            node.lane = i ? cov::ui::NboAomoLane::Centre : cov::ui::NboAomoLane::Left;
            node.quantitative_energy = true;
            node.x = node.label_x = 80.f + 300.f * i;
            node.y = 180.f; node.label_y = 200.f;
            node.width = 100.f; node.height = 10.f;
            node.label_width = 140.f; node.label_height = 20.f;
            view.nodes.push_back(node);
        }
        cov::ui::NboAomoEdge edge;
        edge.id = "fixture-edge"; edge.source_node = 0; edge.target_node = 1;
        edge.coefficient = .9; edge.weight = .81; edge.visible = true;
        view.edges.push_back(edge);
        const auto result = cov::ui::export_nbo_aomo_bundle(view, integration,
            root / "fixture", cov::DiagramExportContent::Images);
        if (!result.svg || !result.png) throw std::runtime_error(result.error);
        std::cout << "Generated all four diagram PNG/SVG pairs in " << root.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
