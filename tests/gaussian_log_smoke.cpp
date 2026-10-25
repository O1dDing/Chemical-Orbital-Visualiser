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

#include "cov/gaussian_log.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_gaussian_log_enrichment.log";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << " Full point group                 C2V     NOp   4\n"
               " Largest Abelian subgroup         C2V     NOp   4\n"
               " Initial guess orbital symmetries:\n"
               "       Occupied  (A1) (A1) (A1) (A1) (A1)\n"
               "       Virtual   (A1) (A1)\n"
               " The electronic state of the initial guess is 1-A1.\n"
               " Orbital symmetries:\n"
               "       Occupied  (A1) (A1) (B2) (A1) (B1)\n"
               "       Virtual   (A1) (B2)\n"
               " The electronic state is 1-A1.\n"
               " S**2 before annihilation     2.1000, after     2.0000\n"
               " SCF Done:  E(UB3LYP) =  -10.0000000000     A.U.\n"
               " The wavefunction has an internal instability.\n"
               " SCF Done:  E(UB3LYP) =  -10.1000000000     A.U.\n"
               " <S^2> before annihilation = 0.7523, after = 0.7500\n"
               " The wavefunction is stable under the perturbations considered.\n";
    }

    cov::Wavefunction wf;
    wf.basis_count = 7;
    wf.orbitals.resize(7);
    for (auto& mo : wf.orbitals) mo.spin = cov::Spin::Alpha;

    const auto result = cov::enrich_from_gaussian_log(wf, path);
    std::error_code ec;
    std::filesystem::remove(path, ec);

    const std::vector<std::string> expected = {
        "A1", "A1", "B2", "A1", "B1", "A1", "B2"
    };
    if (!result.opened || result.symmetry_labels_applied != 7u ||
        !result.point_group_applied || wf.point_group_detected != "C2V" ||
        wf.point_group_used != "C2V" ||
        wf.point_group_provenance != cov::DataProvenance::Producer ||
        !result.spin_squared_applied || !result.scf_diagnostic_applied ||
        !result.stability_diagnostic_applied ||
        wf.scf_convergence != cov::ScfConvergenceStatus::Converged ||
        wf.scf_convergence_provenance != cov::DataProvenance::Producer ||
        wf.stability != cov::WavefunctionStabilityStatus::Stable ||
        wf.stability_provenance != cov::DataProvenance::Producer ||
        wf.spin_squared_provenance != cov::DataProvenance::Producer ||
        std::abs(wf.spin_squared_before_annihilation - 0.7523) > 1.0e-10 ||
        std::abs(wf.spin_squared_after_annihilation - 0.7500) > 1.0e-10) {
        std::cerr << "Gaussian LOG point-group/symmetry metadata regression\n";
        return EXIT_FAILURE;
    }

    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (wf.orbitals[i].symmetry != expected[i] ||
            wf.orbitals[i].symmetry_provenance != cov::DataProvenance::Producer) {
            std::cerr << "Gaussian LOG orbital symmetry ordering regression at " << i << '\n';
            return EXIT_FAILURE;
        }
    }

    // A later explicit SCF failure invalidates S**2 left by an earlier cycle.
    const auto failed_path = std::filesystem::temp_directory_path() /
                             "cov_gaussian_log_failed_scf.log";
    {
        std::ofstream out(failed_path, std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << " S**2 before annihilation     0.9000, after     0.8000\n"
               " SCF Done:  E(UHF) = -1.0 A.U.\n"
               " Convergence failure -- run terminated.\n";
    }
    cov::Wavefunction failed;
    const auto failed_result = cov::enrich_from_gaussian_log(failed, failed_path);
    std::filesystem::remove(failed_path, ec);
    if (!failed_result.scf_diagnostic_applied ||
        failed_result.spin_squared_applied ||
        failed.scf_convergence != cov::ScfConvergenceStatus::Failed ||
        failed.spin_squared_provenance != cov::DataProvenance::Unavailable) {
        std::cerr << "Gaussian LOG failed-SCF/stale-S**2 regression\n";
        return EXIT_FAILURE;
    }

    // Route requests alone are not producer outcomes.
    const auto route_path = std::filesystem::temp_directory_path() /
                            "cov_gaussian_log_route_only.log";
    {
        std::ofstream out(route_path, std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << " #p ub3lyp/def2svp scf=(xqc,tight) stable=opt\n";
    }
    cov::Wavefunction route_only;
    const auto route_result = cov::enrich_from_gaussian_log(route_only, route_path);
    std::filesystem::remove(route_path, ec);
    if (route_result.scf_diagnostic_applied ||
        route_result.stability_diagnostic_applied ||
        route_only.scf_convergence != cov::ScfConvergenceStatus::Unavailable ||
        route_only.stability != cov::WavefunctionStabilityStatus::Unavailable) {
        std::cerr << "Gaussian LOG route-keyword inference regression\n";
        return EXIT_FAILURE;
    }

    // Diagnostics from an earlier SCF cycle must not be attached to a later
    // final state simply because the lines happen to be nearby.
    const auto stale_path = std::filesystem::temp_directory_path() /
                            "cov_gaussian_log_stale_final_state.log";
    {
        std::ofstream out(stale_path, std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << " SCF Done:  E(UHF) = -1.0 A.U.\n"
               " <S^2> before annihilation = 0.8000, after = 0.7500\n"
               " The wavefunction is stable under the perturbations considered.\n"
               " SCF Done:  E(UHF) = -1.1 A.U.\n";
    }
    cov::Wavefunction stale;
    const auto stale_result = cov::enrich_from_gaussian_log(stale, stale_path);
    std::filesystem::remove(stale_path, ec);
    if (!stale_result.scf_diagnostic_applied ||
        stale_result.spin_squared_applied ||
        stale_result.stability_diagnostic_applied ||
        stale.scf_convergence != cov::ScfConvergenceStatus::Converged ||
        stale.spin_squared_provenance != cov::DataProvenance::Unavailable ||
        stale.stability_provenance != cov::DataProvenance::Unavailable) {
        std::cerr << "Gaussian LOG stale-final-state diagnostic regression\n";
        return EXIT_FAILURE;
    }

    // A same-sized orbital block and point group from an earlier Link1 job
    // must not be attached to the final job merely because dimensions match.
    const auto stale_symmetry_path = std::filesystem::temp_directory_path() /
        "cov_gaussian_log_stale_link1_symmetry.log";
    {
        std::ofstream out(stale_symmetry_path, std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << " Entering Link 1 = old-job\n"
               " Full point group C2V NOp 4\n"
               " Largest Abelian subgroup C2V NOp 4\n"
               " Orbital symmetries:\n"
               "   Occupied (A1)\n"
               "   Virtual  (B2)\n"
               " The electronic state is 1-A1.\n"
               " SCF Done: E(RHF) = -1.0 A.U.\n"
               " Entering Link 1 = final-job\n"
               " SCF Done: E(RHF) = -1.1 A.U.\n";
    }
    cov::Wavefunction stale_symmetry;
    stale_symmetry.orbitals.resize(2u);
    const auto stale_symmetry_result=cov::enrich_from_gaussian_log(
        stale_symmetry,stale_symmetry_path);
    std::filesystem::remove(stale_symmetry_path,ec);
    if (stale_symmetry_result.point_group_applied ||
        stale_symmetry_result.symmetry_labels_applied!=0u ||
        stale_symmetry.point_group_provenance!=cov::DataProvenance::Unavailable ||
        std::any_of(stale_symmetry.orbitals.begin(),stale_symmetry.orbitals.end(),
            [](const auto& mo){
                return mo.symmetry_provenance!=cov::DataProvenance::Unavailable;
            })) {
        std::cerr << "Gaussian LOG stale-Link1 symmetry/point-group regression\n";
        return EXIT_FAILURE;
    }

    // The point group is a bounded field, not the next whitespace token.
    // Exercise groups outside the small display catalogue, both linear
    // aliases, invalid/missing fields, and a later empty geometry record.
    const auto field_path=std::filesystem::temp_directory_path()/"cov_gaussian_point_group_field.log";
    for (const std::string value : {"D6H","D12H","C7V","S8","IH","D*H","C*V","KH",""}) {
        {
            std::ofstream out(field_path,std::ios::binary);
            out << " Full point group       " << value << "       NOp 0\n"
                << " Largest Abelian subgroup       " << value << "       NOp 0\n";
        }
        cov::Wavefunction fields;
        const auto enriched=cov::enrich_from_gaussian_log(fields,field_path);
        if (enriched.point_group_applied!=!value.empty() ||
            fields.point_group_detected!=value || fields.point_group_used!=value ||
            fields.point_group_source_records.size()!=2u ||
            fields.point_group_source_records.front().valid!=!value.empty() ||
            fields.point_group_source_records.front().status!=(value.empty()?"missing-input":"available") ||
            fields.point_group_source_records.front().raw_record.find("NOp")==std::string::npos) {
            std::cerr << "Gaussian point-group bounded field regression: " << value << '\n';
            return EXIT_FAILURE;
        }
    }
    for (const std::string value : {"","NOp","banana","C0V","D6H extra"}) {
        {
            std::ofstream out(field_path,std::ios::binary);
            out << " Full point group C2V NOp 4\n"
                << " Full point group " << value << " NOp 0\n"
                << " SCF Done: E(RHF) = -1.0 A.U.\n";
        }
        cov::Wavefunction fields;
        const auto enriched=cov::enrich_from_gaussian_log(fields,field_path);
        if (enriched.point_group_applied || !fields.point_group_detected.empty() ||
            fields.point_group_provenance!=cov::DataProvenance::Unavailable ||
            fields.point_group_source_records.size()!=2u ||
            !fields.point_group_source_records.front().valid ||
            fields.point_group_source_records.back().valid) {
            std::cerr << "Gaussian invalid later field resurrected stale point group\n";
            return EXIT_FAILURE;
        }
    }
    std::filesystem::remove(field_path,ec);

    {
        std::ofstream out(field_path,std::ios::binary);
        out << " Full point group D6H NOp 24\n";
    }
    cov::Wavefunction detected_only;
    (void)cov::enrich_from_gaussian_log(detected_only,field_path);
    std::filesystem::remove(field_path,ec);
    if (detected_only.point_group_detected!="D6H" || !detected_only.point_group_used.empty() ||
        detected_only.point_group_detected_provenance!=cov::DataProvenance::Producer ||
        detected_only.point_group_used_provenance!=cov::DataProvenance::Unavailable) {
        std::cerr << "Gaussian detected group invented an unreported calculation subgroup\n";
        return EXIT_FAILURE;
    }

    // Parallel .chk/.log project trees are a common Gaussian layout. The
    // exact relative subdirectory and stem must be preserved.
    const auto tree_root = std::filesystem::temp_directory_path() /
                           "cov_gaussian_companion_tree";
    const auto fch_path = tree_root / ".chk" / "series" / "frame.fch";
    const auto log_path = tree_root / ".log" / "series" / "frame.log";
    std::filesystem::create_directories(fch_path.parent_path(), ec);
    std::filesystem::create_directories(log_path.parent_path(), ec);
    {
        std::ofstream out(log_path, std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << " SCF Done:  E(RHF) = -1.0 A.U.\n";
    }
    const auto discovered = cov::find_sibling_gaussian_log(fch_path);
    std::filesystem::remove_all(tree_root, ec);
    if (discovered != log_path) {
        std::cerr << "parallel Gaussian .chk/.log companion lookup regression\n";
        return EXIT_FAILURE;
    }

    std::cout << "Gaussian LOG enrichment smoke test passed\n";
    return EXIT_SUCCESS;
}
