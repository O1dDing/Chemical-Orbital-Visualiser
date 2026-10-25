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

#include "cov/coordination_geometry.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include "cov/point_group_catalog.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <span>
#include <string_view>
#include <vector>

namespace {

using Mat3 = std::array<double, 9>;
using Vec3 = cov::CoordinationDirection;

constexpr double kInvSqrt2 = 0.70710678118654752440;
constexpr double kInvSqrt6 = 0.40824829046386301637;

constexpr std::array<Mat3, 5> kDBasis{{
    Mat3{-kInvSqrt6, 0.0, 0.0,
         0.0, -kInvSqrt6, 0.0,
         0.0, 0.0, 2.0 * kInvSqrt6},
    Mat3{0.0, 0.0, kInvSqrt2,
         0.0, 0.0, 0.0,
         kInvSqrt2, 0.0, 0.0},
    Mat3{0.0, 0.0, 0.0,
         0.0, 0.0, kInvSqrt2,
         0.0, kInvSqrt2, 0.0},
    Mat3{kInvSqrt2, 0.0, 0.0,
         0.0, -kInvSqrt2, 0.0,
         0.0, 0.0, 0.0},
    Mat3{0.0, kInvSqrt2, 0.0,
         kInvSqrt2, 0.0, 0.0,
         0.0, 0.0, 0.0},
}};

constexpr Mat3 kC2z{
    -1.0, 0.0, 0.0,
    0.0, -1.0, 0.0,
    0.0, 0.0, 1.0,
};
constexpr Mat3 kSigmaXz{
    1.0, 0.0, 0.0,
    0.0, -1.0, 0.0,
    0.0, 0.0, 1.0,
};
constexpr Mat3 kSigmaYz{
    -1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
};
constexpr Mat3 kC3xyz{
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0, 0.0, 0.0,
};
constexpr Mat3 kS4z{
    0.0, -1.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, -1.0,
};
constexpr Mat3 kC2x{
    1.0, 0.0, 0.0,
    0.0, -1.0, 0.0,
    0.0, 0.0, -1.0,
};
constexpr Mat3 kSigmaDiagonal{
    0.0, 1.0, 0.0,
    1.0, 0.0, 0.0,
    0.0, 0.0, 1.0,
};

Mat3 multiply(const Mat3& left, const Mat3& right) {
    Mat3 result{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            for (std::size_t k = 0; k < 3; ++k) {
                result[3 * row + column] +=
                    left[3 * row + k] * right[3 * k + column];
            }
        }
    }
    return result;
}

Mat3 transpose(const Mat3& value) {
    return {
        value[0], value[3], value[6],
        value[1], value[4], value[7],
        value[2], value[5], value[8],
    };
}

Mat3 rotation_xyz(const double x, const double y, const double z) {
    const double cx = std::cos(x);
    const double sx = std::sin(x);
    const double cy = std::cos(y);
    const double sy = std::sin(y);
    const double cz = std::cos(z);
    const double sz = std::sin(z);
    const Mat3 rx{1.0, 0.0, 0.0, 0.0, cx, -sx, 0.0, sx, cx};
    const Mat3 ry{cy, 0.0, sy, 0.0, 1.0, 0.0, -sy, 0.0, cy};
    const Mat3 rz{cz, -sz, 0.0, sz, cz, 0.0, 0.0, 0.0, 1.0};
    return multiply(rz, multiply(ry, rx));
}

Vec3 transform(const Mat3& matrix, const Vec3& value) {
    return {
        matrix[0] * value.x + matrix[1] * value.y + matrix[2] * value.z,
        matrix[3] * value.x + matrix[4] * value.y + matrix[5] * value.z,
        matrix[6] * value.x + matrix[7] * value.y + matrix[8] * value.z,
    };
}

double squared_distance(const Vec3& left, const Vec3& right) {
    const double x = left.x - right.x;
    const double y = left.y - right.y;
    const double z = left.z - right.z;
    return x * x + y * y + z * z;
}

bool invariant_under(const std::vector<Vec3>& directions,
                     const Mat3& operation,
                     const double tolerance = 2.0e-5) {
    std::vector<bool> used(directions.size(), false);
    for (const auto& direction : directions) {
        const auto moved = transform(operation, direction);
        std::size_t best = directions.size();
        double best_distance = tolerance * tolerance;
        for (std::size_t i = 0; i < directions.size(); ++i) {
            if (used[i]) continue;
            const double distance = squared_distance(moved, directions[i]);
            if (distance <= best_distance) {
                best = i;
                best_distance = distance;
            }
        }
        if (best == directions.size()) return false;
        used[best] = true;
    }
    return true;
}

double operation_residual(const std::vector<Vec3>& directions,
                          const Mat3& operation) {
    double worst_squared = 0.0;
    for (const auto& direction : directions) {
        const auto moved = transform(operation, direction);
        double best_squared = 4.0;
        for (const auto& candidate : directions) {
            best_squared = std::min(
                best_squared, squared_distance(moved, candidate));
        }
        worst_squared = std::max(worst_squared, best_squared);
    }
    return std::sqrt(worst_squared);
}

bool require_operations(const std::string_view machine_id,
                        const std::span<const Mat3> operations,
                        const double tolerance = 2.0e-5) {
    const auto* geometry = cov::coordination_geometry_descriptor(machine_id);
    if (!geometry) {
        std::cerr << machine_id << ": missing descriptor\n";
        return false;
    }
    bool valid = true;
    for (std::size_t i = 0; i < operations.size(); ++i) {
        const auto& operation = operations[i];
        if (!invariant_under(
                geometry->reference_directions, operation, tolerance)) {
            std::cerr << machine_id
                      << ": reference is not invariant under canonical operation #"
                      << (i + 1u) << " (max unit-direction residual "
                      << operation_residual(
                             geometry->reference_directions, operation)
                      << ")\n";
            valid = false;
        }
    }
    return valid;
}

bool validate_catalogue_point_group_contract() {
    for (const auto& geometry : cov::coordination_geometry_catalog()) {
        const auto* group = cov::find_point_group(geometry.point_group);
        const auto decomposition = cov::decompose_metal_spd(geometry.point_group);
        if (!group || group->canonical_name != geometry.point_group ||
            !decomposition ||
            cov::irrep_total_dimension(decomposition->s) != 1u ||
            cov::irrep_total_dimension(decomposition->p) != 3u ||
            cov::irrep_total_dimension(decomposition->d) != 5u ||
            decomposition->total_dimension() != 9u) {
            std::cerr << geometry.machine_id
                      << ": geometry/point-group s-p-d contract failed\n";
            return false;
        }
    }
    return true;
}

bool validate_explicit_canonical_axes() {
    const auto* linear = cov::coordination_geometry_descriptor("L-2");
    if (!linear || linear->reference_directions.size() != 2u) return false;
    for (const auto& direction : linear->reference_directions) {
        if (std::abs(direction.x) > 1.0e-10 ||
            std::abs(direction.y) > 1.0e-10 ||
            std::abs(std::abs(direction.z) - 1.0) > 1.0e-10) {
            std::cerr << "L-2: molecular axis is not canonical z\n";
            return false;
        }
    }
    if (linear->reference_directions[0].z *
            linear->reference_directions[1].z > -0.9999999999) {
        std::cerr << "L-2: endpoints are not antipodal\n";
        return false;
    }

    constexpr std::array<Mat3, 3> c2v_operations{
        kC2z, kSigmaXz, kSigmaYz};
    bool valid = true;
    valid &= require_operations("A-2", c2v_operations);
    valid &= require_operations("TS-3", c2v_operations);
    valid &= require_operations("SS-4", c2v_operations);
    // The published JSPC-10 decimal coordinates are only approximately C2v:
    // rotation-invariant dot products for nominally reflected pairs differ by
    // up to 3.12e-5.  The fitted frame is nevertheless canonical to that
    // source precision, so this single reference receives a 4e-5 allowance.
    valid &= require_operations("SPC-10", c2v_operations, 4.0e-5);

    constexpr std::array<Mat3, 2> tetrahedral_operations{kC2z, kC3xyz};
    valid &= require_operations("T-4", tetrahedral_operations);

    constexpr std::array<Mat3, 3> d2d_operations{
        kS4z, kC2x, kSigmaDiagonal};
    valid &= require_operations("TDD-8", d2d_operations);
    return valid;
}

double matrix_dot(const Mat3& left, const Mat3& right) {
    double result = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) result += left[i] * right[i];
    return result;
}

cov::Wavefunction central_spd_wavefunction() {
    cov::Wavefunction wavefunction;
    wavefunction.atoms.push_back({"M", 24, 0.0, 0.0, 0.0});
    wavefunction.shells.push_back({0, 0, 0, 0, 0, 0});
    wavefunction.shells.push_back({0, 0, 0, 1, 1, 1});
    wavefunction.shells.push_back({0, 0, 0, 4, 2, 1});
    wavefunction.basis_count = 9;
    // This fixture uses orthonormal, single-centre pure angular functions.
    // Supply their metric explicitly; absent overlap data is not identity S.
    wavefunction.ao_overlap.assign(81u, 0.0);
    for (std::size_t i = 0; i < 9u; ++i)
        wavefunction.ao_overlap[i * 9u + i] = 1.0;
    return wavefunction;
}

std::vector<double> rotated_p_coefficients(
    const Mat3& rotation, const std::array<double, 3>& reference) {
    std::array<double, 3> input{};
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            input[row] += rotation[3 * row + column] * reference[column];
        }
    }
    std::vector<double> coefficients(9, 0.0f);
    // Gaussian pure-p order used by the production classifier: pz,-px,-py.
    coefficients[1] = static_cast<float>(input[2]);
    coefficients[2] = static_cast<float>(-input[0]);
    coefficients[3] = static_cast<float>(-input[1]);
    return coefficients;
}

std::vector<double> rotated_d_coefficients(const Mat3& rotation,
                                          const std::size_t component) {
    const Mat3 input = multiply(
        multiply(rotation, kDBasis[component]), transpose(rotation));
    constexpr std::array<double, 5> signs{1.0, -1.0, -1.0, 1.0, 1.0};
    std::vector<double> coefficients(9, 0.0f);
    for (std::size_t i = 0; i < kDBasis.size(); ++i) {
        coefficients[4 + i] = static_cast<float>(
            signs[i] * matrix_dot(input, kDBasis[i]));
    }
    return coefficients;
}

std::size_t add_orbital(cov::Wavefunction& wavefunction,
                        std::vector<double> coefficients) {
    cov::MolecularOrbital orbital;
    orbital.coefficients = std::move(coefficients);
    wavefunction.orbitals.push_back(std::move(orbital));
    return wavefunction.orbitals.size() - 1u;
}

bool expect_label(const cov::Wavefunction& wavefunction,
                  const std::size_t orbital,
                  const std::string_view point_group,
                  const Mat3& matched_frame,
                  const std::string_view expected,
                  const cov::MetalAOShell expected_shell) {
    const std::array<std::size_t, 1> indices{orbital};
    const auto assignment = cov::classify_local_metal_irrep(
        wavefunction, indices, 0, point_group, matched_frame);
    if (!assignment || assignment->label != expected ||
        assignment->shell != expected_shell || assignment->confidence < 0.999) {
        std::cerr << point_group << ": expected " << expected << " from "
                  << (expected_shell == cov::MetalAOShell::P ? 'p' : 'd')
                  << " AO after matcher frame";
        if (assignment) {
            std::cerr << ", got " << assignment->label
                      << " at confidence " << assignment->confidence;
        } else {
            std::cerr << ", got no assignment";
        }
        std::cerr << '\n';
        return false;
    }
    return true;
}

struct CrossFrameCase {
    std::string_view machine_id;
    std::string_view point_group;
    std::string_view p_label;
    std::array<double, 3> p_component;
    std::string_view d_label;
    std::size_t d_component = 0;
};

bool validate_matcher_ao_cross_frames() {
    constexpr std::array<CrossFrameCase, 12> cases{{
        {"L-2", "Dinfh", "Sigma_u+", {0.0, 0.0, 1.0},
         "Sigma_g+", 0},
        {"TP-3", "D3h", "A2''", {0.0, 0.0, 1.0}, "E''", 1},
        {"TPY-3", "C3v", "A1", {0.0, 0.0, 1.0}, "E", 3},
        {"T-4", "Td", "T2", {0.0, 0.0, 1.0}, "E", 0},
        {"SP-4", "D4h", "A2u", {0.0, 0.0, 1.0}, "B2g", 4},
        {"SPY-5", "C4v", "A1", {0.0, 0.0, 1.0}, "B2", 4},
        {"OC-6", "Oh", "T1u", {1.0, 0.0, 0.0}, "T2g", 4},
        {"PBPY-7", "D5h", "A2''", {0.0, 0.0, 1.0}, "E2'", 3},
        {"SAPR-8", "D4d", "E1", {1.0, 0.0, 0.0}, "E3", 1},
        {"TDD-8", "D2d", "B2", {0.0, 0.0, 1.0}, "A1", 0},
        {"PAPR-10", "D5d", "A2u", {0.0, 0.0, 1.0}, "E2g", 3},
        {"SPC-10", "C2v", "A1", {0.0, 0.0, 1.0}, "A2", 4},
    }};
    const std::array<Mat3, 3> global_rotations{
        rotation_xyz(0.37, -0.51, 0.22),
        rotation_xyz(-0.43, 0.29, 0.61),
        rotation_xyz(1.03, -0.17, -0.74),
    };

    for (const auto& item : cases) {
        const auto* geometry =
            cov::coordination_geometry_descriptor(item.machine_id);
        if (!geometry) return false;
        for (const auto& global_rotation : global_rotations) {
            std::vector<Vec3> input;
            input.reserve(geometry->reference_directions.size());
            for (std::size_t i = 0; i < geometry->reference_directions.size();
                 ++i) {
                const auto reference = geometry->reference_directions[
                    geometry->reference_directions.size() - 1u - i];
                auto moved = transform(global_rotation, reference);
                const double radius = 1.7 + 0.03 * static_cast<double>(i);
                moved.x *= radius;
                moved.y *= radius;
                moved.z *= radius;
                input.push_back(moved);
            }
            const auto match = cov::match_coordination_geometry(input);
            if (!match.best || match.best->id != geometry->id ||
                !match.accepted || match.ambiguous ||
                match.best->angular_rms > 2.0e-6) {
                std::cerr << item.machine_id
                          << ": arbitrary-frame matcher did not recover geometry\n";
                return false;
            }

            auto wavefunction = central_spd_wavefunction();
            const auto p = add_orbital(
                wavefunction,
                rotated_p_coefficients(global_rotation, item.p_component));
            const auto d = add_orbital(
                wavefunction,
                rotated_d_coefficients(global_rotation, item.d_component));
            if (!expect_label(
                    wavefunction, p, item.point_group,
                    match.best->rotation_reference_to_input, item.p_label,
                    cov::MetalAOShell::P) ||
                !expect_label(
                    wavefunction, d, item.point_group,
                    match.best->rotation_reference_to_input, item.d_label,
                    cov::MetalAOShell::D)) {
                std::cerr << item.machine_id
                          << ": geometry/AO canonical-frame contract failed\n";
                return false;
            }
        }
    }
    return true;
}

} // namespace

int main() {
    bool valid = true;
    valid &= validate_catalogue_point_group_contract();
    valid &= validate_explicit_canonical_axes();
    valid &= validate_matcher_ao_cross_frames();
    if (!valid) {
        return 1;
    }
    std::cout << "coordination frame contract smoke tests passed\n";
    return 0;
}
