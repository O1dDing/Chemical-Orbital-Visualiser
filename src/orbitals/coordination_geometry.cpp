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

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <utility>

namespace cov {
namespace {

using Direction = CoordinationDirection;
using Matrix3 = std::array<double, 9>;

constexpr Matrix3 kIdentity{
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
};

constexpr double kInvSqrt2 = 0.70710678118654752440;

// Columns are the catalogue's canonical x/y/z axes expressed in the raw
// SHAPE reference frame.  Reference directions are rotated back into the
// shared point_group_catalog convention before matching, so a returned
// reference->input rotation is also the local AO symmetry frame.
constexpr Matrix3 kRawZRotated45{
    kInvSqrt2, -kInvSqrt2, 0.0,
    kInvSqrt2,  kInvSqrt2, 0.0,
    0.0,        0.0,       1.0,
};
constexpr Matrix3 kLinearRawX{
    0.0, 0.0, 1.0,
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
};
constexpr Matrix3 kTShapeRawY{
    1.0,  0.0, 0.0,
    0.0,  0.0, 1.0,
    0.0, -1.0, 0.0,
};
constexpr Matrix3 kSeesawRawDiagonal{
    0.0,       kInvSqrt2, kInvSqrt2,
    0.0,      -kInvSqrt2, kInvSqrt2,
    1.0,       0.0,       0.0,
};
constexpr Matrix3 kSphenocoronaCanonicalToRaw{
     0.9965074934719306,    0.000136096511295039,  0.0835032749778707,
    -0.0835008339282680,   -0.00619391072237890,   0.996488457636739,
     0.000652830432869321, -0.999980808289691,     -0.00616091427655885,
};

double dot(const Direction& left, const Direction& right) noexcept {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

Direction subtract(const Direction& left, const Direction& right) noexcept {
    return {left.x - right.x, left.y - right.y, left.z - right.z};
}

Direction multiply(const Direction& value, const double scale) noexcept {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Direction cross(const Direction& left, const Direction& right) noexcept {
    return {
        left.y * right.z - left.z * right.y,
        left.z * right.x - left.x * right.z,
        left.x * right.y - left.y * right.x,
    };
}

double norm(const Direction& value) noexcept {
    return std::sqrt(dot(value, value));
}

Direction normalised(const Direction& value) noexcept {
    const double length = norm(value);
    if (!(length > 1.0e-12) || !std::isfinite(length)) return {};
    return multiply(value, 1.0 / length);
}

Direction rotate_vector(const Matrix3& matrix, const Direction& value) noexcept {
    return {
        matrix[0] * value.x + matrix[1] * value.y + matrix[2] * value.z,
        matrix[3] * value.x + matrix[4] * value.y + matrix[5] * value.z,
        matrix[6] * value.x + matrix[7] * value.y + matrix[8] * value.z,
    };
}

Matrix3 transpose(const Matrix3& matrix) noexcept {
    return {
        matrix[0], matrix[3], matrix[6],
        matrix[1], matrix[4], matrix[7],
        matrix[2], matrix[5], matrix[8],
    };
}

std::vector<Direction> centred_unit_directions(
    const std::initializer_list<Direction> ligand_positions,
    const Direction centre = {}) {
    std::vector<Direction> result;
    result.reserve(ligand_positions.size());
    for (const auto& position : ligand_positions) {
        result.push_back(normalised(subtract(position, centre)));
    }
    return result;
}

std::vector<Direction> in_symmetry_frame(
    std::vector<Direction> raw_directions,
    const Matrix3& canonical_to_raw) {
    const Matrix3 raw_to_canonical=transpose(canonical_to_raw);
    for (auto& direction:raw_directions) {
        direction=rotate_vector(raw_to_canonical,direction);
    }
    return raw_directions;
}

CoordinationGeometryDescriptor descriptor(
    const GeometryId id,
    const std::string_view machine_id,
    const std::string_view name,
    const std::string_view point_group,
    const std::size_t coordination_number,
    std::vector<Direction> reference_directions,
    const std::string_view reference_note) {
    CoordinationGeometryDescriptor result;
    result.id = id;
    result.machine_id = machine_id;
    result.name = name;
    result.point_group = point_group;
    result.coordination_number = coordination_number;
    result.reference_status = ReferenceGeometryStatus::Available;
    result.reference_note = reference_note;
    result.reference_directions = std::move(reference_directions);
    return result;
}

const std::vector<CoordinationGeometryDescriptor>& catalog_storage() {
    // The available references below are metal-centred unit directions made
    // from the CoSyMlib ideal_structures_center coordinate transcription used
    // by q-shape (MIT, HenriqueCSJ/q-shape).  The source is recorded on every
    // descriptor so replacing a reference is an auditable data change.
    static const std::vector<CoordinationGeometryDescriptor> catalog = [] {
        constexpr std::string_view cosymlib_reference =
            "SHAPE/CoSyMlib ideal_structures_center reference; converted to metal-centred unit directions";
        std::vector<CoordinationGeometryDescriptor> values;
        values.reserve(26);

        values.push_back(descriptor(
            GeometryId::Linear2, "L-2", "Linear", "Dinfh", 2,
            in_symmetry_frame(centred_unit_directions({
                {1.224744871392, 0.0, 0.0},
                {-1.224744871392, 0.0, 0.0},
            }),kLinearRawX), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::Angular2, "A-2", "Angular", "C2v", 2,
            in_symmetry_frame(centred_unit_directions({
                {0.801783725737, 0.801783725737, 0.267261241912},
                {-0.801783725737, -0.801783725737, 0.267261241912},
            }, {0.0, 0.0, -0.534522483825}),kRawZRotated45),
            "SHAPE vT-2 (109.47 degree angular) reference; exposed as COV A-2"));

        values.push_back(descriptor(
            GeometryId::TrigonalPlanar3, "TP-3", "Trigonal planar", "D3h", 3,
            centred_unit_directions({
                {1.154700538379, 0.0, 0.0},
                {-0.577350269190, 1.0, 0.0},
                {-0.577350269190, -1.0, 0.0},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::TrigonalPyramidal3, "TPY-3", "Trigonal pyramidal", "C3v", 3,
            centred_unit_directions({
                {1.137070487230, 0.0, 0.100503781526},
                {-0.568535243615, 0.984731927835, 0.100503781526},
                {-0.568535243615, -0.984731927835, 0.100503781526},
            }, {0.0, 0.0, -0.301511344578}),
            "SHAPE vT-3 reference; exposed as COV TPY-3"));
        values.push_back(descriptor(
            GeometryId::TShaped3, "TS-3", "T-shaped", "C2v", 3,
            in_symmetry_frame(centred_unit_directions({
                {1.206045378311, -0.301511344578, 0.0},
                {0.0, 0.904534033733, 0.0},
                {-1.206045378311, -0.301511344578, 0.0},
            }, {0.0, -0.301511344578, 0.0}),kTShapeRawY),
            "SHAPE mer-vOC-3 reference; exposed as COV TS-3"));

        values.push_back(descriptor(
            GeometryId::Tetrahedral4, "T-4", "Tetrahedral", "Td", 4,
            in_symmetry_frame(centred_unit_directions({
                {0.0, 0.912870929175, -0.645497224368},
                {0.0, -0.912870929175, -0.645497224368},
                {0.912870929175, 0.0, 0.645497224368},
                {-0.912870929175, 0.0, 0.645497224368},
            }),kRawZRotated45), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::SquarePlanar4, "SP-4", "Square planar", "D4h", 4,
            centred_unit_directions({
                {1.118033988750, 0.0, 0.0},
                {0.0, 1.118033988750, 0.0},
                {-1.118033988750, 0.0, 0.0},
                {0.0, -1.118033988750, 0.0},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::Seesaw4, "SS-4", "Seesaw", "C2v", 4,
            in_symmetry_frame(centred_unit_directions({
                {-0.235702260396, -0.235702260396, -1.178511301978},
                {0.942809041582, -0.235702260396, 0.0},
                {-0.235702260396, 0.942809041582, 0.0},
                {-0.235702260396, -0.235702260396, 1.178511301978},
            }, {-0.235702260396, -0.235702260396, 0.0}),
            kSeesawRawDiagonal),
            cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::TrigonalPyramidal4, "vTBPY-4",
            "Trigonal-pyramidal (axially vacant trigonal bipyramid)",
            "C3v", 4,
            centred_unit_directions({
                {0.0, 0.0, -0.917662935482},
                {1.147078669353, 0.0, 0.229415733871},
                {-0.573539334676, 0.993399267799, 0.229415733871},
                {-0.573539334676, -0.993399267799, 0.229415733871},
            }, {0.0, 0.0, 0.229415733871}),cosymlib_reference));

        values.push_back(descriptor(
            GeometryId::TrigonalBipyramidal5, "TBPY-5", "Trigonal bipyramidal", "D3h", 5,
            centred_unit_directions({
                {0.0, 0.0, -1.095445115010},
                {1.095445115010, 0.0, 0.0},
                {-0.547722557505, 0.948683298051, 0.0},
                {-0.547722557505, -0.948683298051, 0.0},
                {0.0, 0.0, 1.095445115010},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::SquarePyramidal5, "SPY-5", "Square pyramidal", "C4v", 5,
            centred_unit_directions({
                {0.0, 0.0, 1.095445115010},
                {1.060660171780, 0.0, -0.273861278753},
                {0.0, 1.060660171780, -0.273861278753},
                {-1.060660171780, 0.0, -0.273861278753},
                {0.0, -1.060660171780, -0.273861278753},
            }), cosymlib_reference));

        values.push_back(descriptor(
            GeometryId::Octahedral6, "OC-6", "Octahedral", "Oh", 6,
            centred_unit_directions({
                {0.0, 0.0, -1.080123449735},
                {1.080123449735, 0.0, 0.0},
                {0.0, 1.080123449735, 0.0},
                {-1.080123449735, 0.0, 0.0},
                {0.0, -1.080123449735, 0.0},
                {0.0, 0.0, 1.080123449735},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::TrigonalPrismatic6, "TPR-6", "Trigonal prismatic", "D3h", 6,
            centred_unit_directions({
                {0.816496580928, 0.0, -0.707106781187},
                {-0.408248290464, 0.707106781187, -0.707106781187},
                {-0.408248290464, -0.707106781187, -0.707106781187},
                {0.816496580928, 0.0, 0.707106781187},
                {-0.408248290464, 0.707106781187, 0.707106781187},
                {-0.408248290464, -0.707106781187, 0.707106781187},
            }), cosymlib_reference));

        values.push_back(descriptor(
            GeometryId::PentagonalBipyramidal7, "PBPY-7", "Pentagonal bipyramidal", "D5h", 7,
            centred_unit_directions({
                {0.0, 0.0, -1.069044967650},
                {1.069044967650, 0.0, 0.0},
                {0.330353062755, 1.016722182696, 0.0},
                {-0.864875546580, 0.628368866022, 0.0},
                {-0.864875546580, -0.628368866022, 0.0},
                {0.330353062755, -1.016722182696, 0.0},
                {0.0, 0.0, 1.069044967650},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::CappedOctahedral7, "COC-7", "Capped octahedral", "C3v", 7,
            centred_unit_directions({
                {0.0, 0.0, 1.128906708829},
                {0.0, -1.046937018035, 0.283078548570},
                {0.906674032650, 0.523468509017, 0.283078548570},
                {-0.906674032650, 0.523468509017, 0.283078548570},
                {0.672964536915, -0.388536257092, -0.678734552207},
                {-0.672964536915, -0.388536257092, -0.678734552207},
                {0.0, 0.777072514184, -0.678734552207},
            }, {0.0, 0.0, 0.058061302083}), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::CappedTrigonalPrismatic7, "CTPR-7", "Capped trigonal prismatic", "C2v", 7,
            centred_unit_directions({
                {0.0, 0.0, 1.020027096827},
                {0.735247575071, 0.735247575071, 0.203750780644},
                {-0.735247575071, 0.735247575071, 0.203750780644},
                {0.735247575071, -0.735247575071, 0.203750780644},
                {-0.735247575071, -0.735247575071, 0.203750780644},
                {0.660960557032, 0.0, -0.892328424325},
                {-0.660960557032, 0.0, -0.892328424325},
            }, {0.0, 0.0, -0.050373370753}), cosymlib_reference));

        values.push_back(descriptor(
            GeometryId::SquareAntiprismatic8, "SAPR-8", "Square antiprismatic", "D4d", 8,
            centred_unit_directions({
                {0.644649377827, 0.644649377827, -0.542083350910},
                {-0.644649377827, 0.644649377827, -0.542083350910},
                {-0.644649377827, -0.644649377827, -0.542083350910},
                {0.644649377827, -0.644649377827, -0.542083350910},
                {0.911671893098, 0.0, 0.542083350910},
                {0.0, 0.911671893098, 0.542083350910},
                {-0.911671893098, 0.0, 0.542083350910},
                {0.0, -0.911671893098, 0.542083350910},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::TriangularDodecahedral8, "TDD-8", "Triangular dodecahedral", "D2d", 8,
            in_symmetry_frame(centred_unit_directions({
                {-0.636106245143, 0.0, 0.848768388024},
                {-0.000000009579, -0.993210924257, 0.372146720241},
                {0.636106254722, 0.0, 0.848768388024},
                {-0.000000009579, 0.993210924257, 0.372146720241},
                {-0.993210876363, 0.0, -0.372146742591},
                {-0.000000009579, -0.636106206828, -0.848768374454},
                {0.993210914678, 0.0, -0.372146742591},
                {-0.000000009579, 0.636106206828, -0.848768374454},
            }, {-0.000000009579, 0.0, 0.000000017561}),
            kRawZRotated45), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::BicappedTrigonalPrismatic8, "BTPR-8", "Bicapped trigonal prismatic", "C2v", 8,
            centred_unit_directions({
                {0.699237877649, 0.0, 0.688732178156},
                {-0.699237877649, 0.0, 0.688732178156},
                {0.699237877649, 0.699237877649, -0.522383347216},
                {-0.699237877649, 0.699237877649, -0.522383347216},
                {0.699237877649, -0.699237877649, -0.522383347216},
                {-0.699237877649, -0.699237877649, -0.522383347216},
                {0.0, 0.925004726938, 0.415373590668},
                {0.0, -0.925004726938, 0.415373590668},
            }, {0.0, 0.0, -0.118678148784}), cosymlib_reference));

        values.push_back(descriptor(
            GeometryId::CappedSquareAntiprismatic9, "CSAPR-9", "Capped square antiprismatic", "C4v", 9,
            centred_unit_directions({
                {0.0, 0.0, 1.053083142672},
                {0.982653581851, 0.0, 0.380440156580},
                {0.0, 0.982653581851, 0.380440156580},
                {-0.982653581851, 0.0, 0.380440156580},
                {0.0, -0.982653581851, 0.380440156580},
                {0.590919690170, 0.590919690170, -0.643458455172},
                {-0.590919690170, 0.590919690170, -0.643458455172},
                {-0.590919690170, -0.590919690170, -0.643458455172},
                {0.590919690170, -0.590919690170, -0.643458455172},
            }, {0.0, 0.0, -0.001009948303}), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::TricappedTrigonalPrismatic9, "TCTPR-9", "Tricapped trigonal prismatic", "D3h", 9,
            centred_unit_directions({
                {0.702728368926, 0.0, 0.785674201318},
                {-0.351364184463, 0.608580619450, 0.785674201318},
                {-0.351364184463, -0.608580619450, 0.785674201318},
                {0.702728368926, 0.0, -0.785674201318},
                {-0.351364184463, 0.608580619450, -0.785674201318},
                {-0.351364184463, -0.608580619450, -0.785674201318},
                {-1.054092553389, 0.0, 0.0},
                {0.527046276695, 0.912870929175, 0.0},
                {0.527046276695, -0.912870929175, 0.0},
            }), cosymlib_reference));

        values.push_back(descriptor(
            GeometryId::PentagonalPrismatic10, "PPR-10", "Pentagonal prism", "D5h", 10,
            centred_unit_directions({
                {0.904182012090, 0.0, -0.531464852095},
                {0.279407607744, 0.859928194515, -0.531464852095},
                {-0.731498613789, 0.531464852095, -0.531464852095},
                {-0.731498613789, -0.531464852095, -0.531464852095},
                {0.279407607744, -0.859928194515, -0.531464852095},
                {0.904182012090, 0.0, 0.531464852095},
                {0.279407607744, 0.859928194515, 0.531464852095},
                {-0.731498613789, 0.531464852095, 0.531464852095},
                {-0.731498613789, -0.531464852095, 0.531464852095},
                {0.279407607744, -0.859928194515, 0.531464852095},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::PentagonalAntiprismatic10, "PAPR-10",
            "Pentagonal antiprism", "D5d", 10,
            centred_unit_directions({
                {0.758925212076, 0.551391442149, -0.469041575982},
                {-0.289883636094, 0.892170094503, -0.469041575982},
                {-0.938083151965, 0.0, -0.469041575982},
                {-0.289883636094, -0.892170094503, -0.469041575982},
                {0.758925212076, -0.551391442149, -0.469041575982},
                {0.938083151965, 0.0, 0.469041575982},
                {0.289883636094, 0.892170094503, 0.469041575982},
                {-0.758925212076, 0.551391442149, 0.469041575982},
                {-0.758925212076, -0.551391442149, 0.469041575982},
                {0.289883636094, -0.892170094503, 0.469041575982},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::BicappedSquareAntiprismatic10, "BCSAPR-10", "Bicapped square antiprismatic", "D4d", 10,
            centred_unit_directions({
                {0.831394933130, 0.0, 0.494350384928},
                {0.587884995060, 0.587884995060, -0.494350384928},
                {0.0, 0.831394933130, 0.494350384928},
                {-0.587884995060, 0.587884995060, -0.494350384928},
                {-0.831394933130, 0.0, 0.494350384928},
                {-0.587884995060, -0.587884995060, -0.494350384928},
                {0.0, -0.831394933130, 0.494350384928},
                {0.587884995060, -0.587884995060, -0.494350384928},
                {0.0, 0.0, 1.325745318058},
                {0.0, 0.0, -1.325745318058},
            }), cosymlib_reference));
        values.push_back(descriptor(
            GeometryId::Sphenocorona10, "SPC-10", "Sphenocorona", "C2v", 10,
            in_symmetry_frame(centred_unit_directions({
                {-1.001871872522, -0.083830395389, -0.581156124487},
                {-1.002034699262, -0.076631127394, 0.581868755803},
                {-0.516334164993, 0.802168048137, -0.005028597236},
                {0.028693454961, 0.335227480386, -0.920231179588},
                {-0.064315504895, -0.772040872012, -0.576759802513},
                {-0.064478331635, -0.764829973537, 0.586265077777},
                {0.028437584370, 0.346602091208, 0.916012486785},
                {0.642643307766, 0.705053528342, -0.004284246426},
                {0.974705182575, -0.249460081184, -0.579853510569},
                {0.974553986317, -0.242260813190, 0.583171369721},
            }, {0.000001057316, 0.000002114633, -0.000004229266}),
            kSphenocoronaCanonicalToRaw),
            "Q-Shape/CoSyMlib SHAPE JSPC-10 reference (HenriqueCSJ/q-shape, MIT)"));
        values.push_back(descriptor(
            GeometryId::Tetradecahedral10, "TD-10", "Tetradecahedral (2:6:2)", "C2v", 10,
            centred_unit_directions({
                {-0.524413653847, 0.908284448462, 0.0},
                {0.524413653847, 0.908284448462, 0.0},
                {-1.048827307693, 0.0, 0.0},
                {1.048827307693, 0.0, 0.0},
                {-0.524413653847, -0.908284448462, 0.0},
                {0.524413653847, -0.908284448462, 0.0},
                {-0.524413653847, 0.0, 0.908284448462},
                {0.524413653847, 0.0, 0.908284448462},
                {0.0, 0.524413653847, -0.908284448462},
                {0.0, -0.524413653847, -0.908284448462},
            }),
            "Q-Shape/CoSyMlib SHAPE TD-10 reference (HenriqueCSJ/q-shape, MIT); point group follows SHAPE C2v table"));
        return values;
    }();
    return catalog;
}

std::vector<std::size_t> minimum_assignment(
    const std::vector<Direction>& input,
    const std::vector<Direction>& reference,
    const Matrix3& input_to_reference) {
    const std::size_t count = input.size();
    std::vector<double> u(count + 1, 0.0);
    std::vector<double> v(count + 1, 0.0);
    std::vector<std::size_t> p(count + 1, 0);
    std::vector<std::size_t> way(count + 1, 0);

    const auto cost = [&](const std::size_t input_index,
                          const std::size_t reference_index) {
        const auto rotated = rotate_vector(input_to_reference, input[input_index]);
        return std::max(0.0, 2.0 - 2.0 * dot(rotated, reference[reference_index]));
    };

    for (std::size_t i = 1; i <= count; ++i) {
        p[0] = i;
        std::size_t column = 0;
        std::vector<double> minimum(count + 1,
                                    std::numeric_limits<double>::infinity());
        std::vector<bool> used(count + 1, false);
        do {
            used[column] = true;
            const std::size_t row = p[column];
            double delta = std::numeric_limits<double>::infinity();
            std::size_t next_column = 0;
            for (std::size_t j = 1; j <= count; ++j) {
                if (used[j]) continue;
                const double reduced = cost(row - 1, j - 1) - u[row] - v[j];
                if (reduced < minimum[j]) {
                    minimum[j] = reduced;
                    way[j] = column;
                }
                if (minimum[j] < delta) {
                    delta = minimum[j];
                    next_column = j;
                }
            }
            for (std::size_t j = 0; j <= count; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minimum[j] -= delta;
                }
            }
            column = next_column;
        } while (p[column] != 0);

        do {
            const std::size_t previous = way[column];
            p[column] = p[previous];
            column = previous;
        } while (column != 0);
    }

    std::vector<std::size_t> assignment(count, 0);
    for (std::size_t column = 1; column <= count; ++column) {
        assignment[p[column] - 1] = column - 1;
    }
    return assignment;
}

Matrix3 quaternion_rotation(const std::array<double, 4>& q) noexcept {
    const double w = q[0];
    const double x = q[1];
    const double y = q[2];
    const double z = q[3];
    return {
        1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z),
        2.0 * (x * z + w * y),
        2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z),
        2.0 * (y * z - w * x),
        2.0 * (x * z - w * y), 2.0 * (y * z + w * x),
        1.0 - 2.0 * (x * x + y * y),
    };
}

std::array<double, 4> largest_eigenvector_4x4(
    std::array<std::array<double, 4>, 4> matrix) {
    std::array<std::array<double, 4>, 4> eigenvectors{};
    for (std::size_t i = 0; i < 4; ++i) eigenvectors[i][i] = 1.0;

    for (std::size_t iteration = 0; iteration < 64; ++iteration) {
        std::size_t p = 0;
        std::size_t q = 1;
        double largest = std::abs(matrix[p][q]);
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) {
                const double magnitude = std::abs(matrix[i][j]);
                if (magnitude > largest) {
                    largest = magnitude;
                    p = i;
                    q = j;
                }
            }
        }
        if (largest < 1.0e-13) break;

        const double angle = 0.5 * std::atan2(
            2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        const double app = matrix[p][p];
        const double aqq = matrix[q][q];
        const double apq = matrix[p][q];
        matrix[p][p] = cosine * cosine * app -
                       2.0 * sine * cosine * apq + sine * sine * aqq;
        matrix[q][q] = sine * sine * app +
                       2.0 * sine * cosine * apq + cosine * cosine * aqq;
        matrix[p][q] = 0.0;
        matrix[q][p] = 0.0;
        for (std::size_t k = 0; k < 4; ++k) {
            if (k == p || k == q) continue;
            const double akp = matrix[k][p];
            const double akq = matrix[k][q];
            matrix[k][p] = cosine * akp - sine * akq;
            matrix[p][k] = matrix[k][p];
            matrix[k][q] = sine * akp + cosine * akq;
            matrix[q][k] = matrix[k][q];
        }
        for (std::size_t k = 0; k < 4; ++k) {
            const double vkp = eigenvectors[k][p];
            const double vkq = eigenvectors[k][q];
            eigenvectors[k][p] = cosine * vkp - sine * vkq;
            eigenvectors[k][q] = sine * vkp + cosine * vkq;
        }
    }

    std::size_t largest_index = 0;
    for (std::size_t i = 1; i < 4; ++i) {
        if (matrix[i][i] > matrix[largest_index][largest_index]) {
            largest_index = i;
        }
    }
    std::array<double, 4> result{};
    double length_squared = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
        result[i] = eigenvectors[i][largest_index];
        length_squared += result[i] * result[i];
    }
    const double inverse_length = 1.0 / std::sqrt(std::max(1.0e-30, length_squared));
    for (auto& value : result) value *= inverse_length;
    return result;
}

Matrix3 optimal_input_to_reference_rotation(
    const std::vector<Direction>& input,
    const std::vector<Direction>& reference,
    const std::vector<std::size_t>& assignment) {
    double sxx = 0.0, sxy = 0.0, sxz = 0.0;
    double syx = 0.0, syy = 0.0, syz = 0.0;
    double szx = 0.0, szy = 0.0, szz = 0.0;
    for (std::size_t i = 0; i < input.size(); ++i) {
        const auto& source = input[i];
        const auto& target = reference[assignment[i]];
        sxx += source.x * target.x;
        sxy += source.x * target.y;
        sxz += source.x * target.z;
        syx += source.y * target.x;
        syy += source.y * target.y;
        syz += source.y * target.z;
        szx += source.z * target.x;
        szy += source.z * target.y;
        szz += source.z * target.z;
    }

    std::array<std::array<double, 4>, 4> davenport{};
    davenport[0] = {sxx + syy + szz, syz - szy, szx - sxz, sxy - syx};
    davenport[1] = {syz - szy, sxx - syy - szz, sxy + syx, szx + sxz};
    davenport[2] = {szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy};
    davenport[3] = {sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz};
    return quaternion_rotation(largest_eigenvector_4x4(davenport));
}

std::optional<Matrix3> pair_frame_rotation(
    const Direction& source_first,
    const Direction& source_second,
    const Direction& target_first,
    const Direction& target_second) {
    const Direction source_axis = normalised(source_first);
    const Direction target_axis = normalised(target_first);
    const Direction source_plane = normalised(subtract(
        source_second, multiply(source_axis, dot(source_axis, source_second))));
    const Direction target_plane = normalised(subtract(
        target_second, multiply(target_axis, dot(target_axis, target_second))));
    if (norm(source_plane) < 0.5 || norm(target_plane) < 0.5) return std::nullopt;
    const Direction source_normal = cross(source_axis, source_plane);
    const Direction target_normal = cross(target_axis, target_plane);
    Matrix3 rotation{};
    const std::array<Direction, 3> source_basis{
        source_axis, source_plane, source_normal};
    const std::array<Direction, 3> target_basis{
        target_axis, target_plane, target_normal};
    for (std::size_t row = 0; row < 3; ++row) {
        const std::array<double, 3> target_components{
            row == 0 ? target_basis[0].x : row == 1 ? target_basis[0].y : target_basis[0].z,
            row == 0 ? target_basis[1].x : row == 1 ? target_basis[1].y : target_basis[1].z,
            row == 0 ? target_basis[2].x : row == 1 ? target_basis[2].y : target_basis[2].z,
        };
        for (std::size_t column = 0; column < 3; ++column) {
            const std::array<double, 3> source_components{
                column == 0 ? source_basis[0].x : column == 1 ? source_basis[0].y : source_basis[0].z,
                column == 0 ? source_basis[1].x : column == 1 ? source_basis[1].y : source_basis[1].z,
                column == 0 ? source_basis[2].x : column == 1 ? source_basis[2].y : source_basis[2].z,
            };
            rotation[row * 3 + column] =
                target_components[0] * source_components[0] +
                target_components[1] * source_components[1] +
                target_components[2] * source_components[2];
        }
    }
    return rotation;
}

struct Alignment {
    Matrix3 input_to_reference = kIdentity;
    std::vector<std::size_t> assignment;
    double angular_rms = std::numeric_limits<double>::infinity();
};

Alignment refine_alignment(
    const std::vector<Direction>& input,
    const std::vector<Direction>& reference,
    Matrix3 rotation,
    const std::size_t max_iterations) {
    std::vector<std::size_t> previous;
    for (std::size_t iteration = 0;
         iteration < std::max<std::size_t>(1, max_iterations); ++iteration) {
        auto assignment = minimum_assignment(input, reference, rotation);
        const Matrix3 fitted = optimal_input_to_reference_rotation(
            input, reference, assignment);
        const bool stable = assignment == previous;
        previous = std::move(assignment);
        rotation = fitted;
        if (stable) break;
    }
    auto assignment = minimum_assignment(input, reference, rotation);
    rotation = optimal_input_to_reference_rotation(input, reference, assignment);
    assignment = minimum_assignment(input, reference, rotation);

    double squared = 0.0;
    for (std::size_t i = 0; i < input.size(); ++i) {
        const double cosine = std::clamp(
            dot(rotate_vector(rotation, input[i]), reference[assignment[i]]),
            -1.0, 1.0);
        const double angle = std::acos(cosine);
        squared += angle * angle;
    }
    return {
        rotation,
        std::move(assignment),
        std::sqrt(squared / static_cast<double>(input.size())),
    };
}

Alignment best_alignment(
    const std::vector<Direction>& input,
    const std::vector<Direction>& reference,
    const std::size_t max_iterations) {
    Alignment best = refine_alignment(input, reference, kIdentity, max_iterations);
    const std::size_t count = input.size();
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = 0; j < count; ++j) {
            if (i == j || std::abs(dot(input[i], input[j])) > 0.995) continue;
            std::array<std::pair<double, std::pair<std::size_t, std::size_t>>, 3>
                closest{{
                    {std::numeric_limits<double>::infinity(), {0, 0}},
                    {std::numeric_limits<double>::infinity(), {0, 0}},
                    {std::numeric_limits<double>::infinity(), {0, 0}},
                }};
            const double input_dot = dot(input[i], input[j]);
            for (std::size_t a = 0; a < count; ++a) {
                for (std::size_t b = 0; b < count; ++b) {
                    if (a == b || std::abs(dot(reference[a], reference[b])) > 0.995) {
                        continue;
                    }
                    const auto candidate = std::make_pair(
                        std::abs(input_dot - dot(reference[a], reference[b])),
                        std::make_pair(a, b));
                    if (candidate < closest[2]) {
                        closest[2] = candidate;
                        if (closest[2] < closest[1]) std::swap(closest[2], closest[1]);
                        if (closest[1] < closest[0]) std::swap(closest[1], closest[0]);
                    }
                }
            }
            for (const auto& seed : closest) {
                if (!std::isfinite(seed.first)) continue;
                const auto rotation = pair_frame_rotation(
                    input[i], input[j], reference[seed.second.first],
                    reference[seed.second.second]);
                if (!rotation) continue;
                auto alignment = refine_alignment(
                    input, reference, *rotation, max_iterations);
                if (alignment.angular_rms + 1.0e-12 < best.angular_rms) {
                    best = std::move(alignment);
                }
            }
        }
    }
    return best;
}

double input_radial_cv(const std::vector<Direction>& input_vectors) {
    std::vector<double> radii;
    radii.reserve(input_vectors.size());
    for (const auto& vector : input_vectors) radii.push_back(norm(vector));
    const double mean = std::accumulate(radii.begin(), radii.end(), 0.0) /
                        static_cast<double>(radii.size());
    if (!(mean > 1.0e-12)) return 0.0;
    double squared = 0.0;
    for (const double radius : radii) {
        const double delta = radius - mean;
        squared += delta * delta;
    }
    return std::sqrt(squared / static_cast<double>(radii.size())) / mean;
}

GeometryCandidate evaluate_candidate(
    const std::vector<Direction>& input,
    const double radial_cv,
    const CoordinationGeometryDescriptor& geometry,
    const GeometryMatchOptions& options) {
    const auto alignment = best_alignment(
        input, geometry.reference_directions, options.max_iterations);
    GeometryCandidate candidate;
    candidate.id = geometry.id;
    candidate.angular_rms = alignment.angular_rms;
    candidate.shape_measure = 100.0 * alignment.angular_rms * alignment.angular_rms;
    candidate.radial_cv = radial_cv;
    candidate.rotation_reference_to_input = transpose(alignment.input_to_reference);
    candidate.assignment = alignment.assignment;
    return candidate;
}

double direction_dot(const CoordinationContact& left,
                     const CoordinationContact& right) noexcept {
    return left.direction.x * right.direction.x +
           left.direction.y * right.direction.y +
           left.direction.z * right.direction.z;
}

CoordinationContact make_contact(
    const Wavefunction& wavefunction,
    const std::size_t centre_atom,
    const std::size_t other_atom,
    const double mayer) {
    const auto& centre = wavefunction.atoms[centre_atom];
    const auto& other = wavefunction.atoms[other_atom];
    const Direction displacement{
        other.x - centre.x,
        other.y - centre.y,
        other.z - centre.z,
    };
    CoordinationContact result;
    result.atom_index = other_atom;
    result.distance = norm(displacement);
    result.mayer_bond_order = std::abs(mayer);
    result.direction = normalised(displacement);
    return result;
}

} // namespace

const std::vector<CoordinationGeometryDescriptor>&
coordination_geometry_catalog() {
    return catalog_storage();
}

const CoordinationGeometryDescriptor*
coordination_geometry_descriptor(const GeometryId id) noexcept {
    const auto& catalog = catalog_storage();
    const auto found = std::find_if(catalog.begin(), catalog.end(),
                                    [&](const auto& item) { return item.id == id; });
    return found == catalog.end() ? nullptr : &*found;
}

const CoordinationGeometryDescriptor*
coordination_geometry_descriptor(const std::string_view machine_id) noexcept {
    const auto& catalog = catalog_storage();
    const auto found = std::find_if(
        catalog.begin(), catalog.end(),
        [&](const auto& item) { return item.machine_id == machine_id; });
    return found == catalog.end() ? nullptr : &*found;
}

std::vector<const CoordinationGeometryDescriptor*>
coordination_geometries_for_cn(const std::size_t coordination_number) {
    std::vector<const CoordinationGeometryDescriptor*> result;
    for (const auto& geometry : catalog_storage()) {
        if (geometry.coordination_number == coordination_number) {
            result.push_back(&geometry);
        }
    }
    return result;
}

GeometryMatch match_coordination_geometry(
    const std::vector<CoordinationDirection>& input_vectors,
    const GeometryMatchOptions& options) {
    GeometryMatch result;
    result.coordination_number = input_vectors.size();
    if (input_vectors.size() < 2 || input_vectors.size() > 10) return result;

    std::vector<Direction> input;
    input.reserve(input_vectors.size());
    for (const auto& vector : input_vectors) {
        if (!std::isfinite(vector.x) || !std::isfinite(vector.y) ||
            !std::isfinite(vector.z) || norm(vector) <= 1.0e-12) {
            return result;
        }
        input.push_back(normalised(vector));
    }
    const double radial_cv = input_radial_cv(input_vectors);

    std::vector<GeometryCandidate> candidates;
    for (const auto* geometry : coordination_geometries_for_cn(input.size())) {
        if (!geometry->matchable()) {
            result.unavailable_candidates.push_back(geometry->id);
            continue;
        }
        candidates.push_back(evaluate_candidate(input, radial_cv, *geometry, options));
    }
    std::sort(candidates.begin(), candidates.end(), [](const auto& left,
                                                       const auto& right) {
        if (left.angular_rms != right.angular_rms) {
            return left.angular_rms < right.angular_rms;
        }
        return static_cast<int>(left.id) < static_cast<int>(right.id);
    });
    if (candidates.empty()) return result;

    result.best = candidates.front();
    if (candidates.size() > 1) result.runner_up = candidates[1];
    result.accepted = result.best->angular_rms <=
                      std::max(0.0, options.max_angular_rms);
    result.ambiguous = result.runner_up.has_value() &&
        result.runner_up->angular_rms - result.best->angular_rms <=
            std::max(0.0, options.ambiguity_margin);
    return result;
}

GeometryMatch analyse_coordination_shell(
    const CoordinationShell& shell,
    const GeometryMatchOptions& options) {
    std::vector<Direction> displacements;
    displacements.reserve(shell.contacts.size());
    for (const auto& contact : shell.contacts) {
        displacements.push_back(multiply(contact.direction, contact.distance));
    }
    return match_coordination_geometry(displacements, options);
}

CoordinationShell extract_coordination_shell(
    const Wavefunction& wavefunction,
    const std::size_t centre_atom,
    const CoordinationShellOptions& options) {
    CoordinationShell result;
    result.centre_atom = centre_atom;
    if (centre_atom >= wavefunction.atoms.size()) return result;

    std::vector<std::pair<std::size_t, double>> bonded;
    for (const auto& bond : wavefunction.bond_orders) {
        std::size_t other = wavefunction.atoms.size();
        if (bond.atom_a == centre_atom) other = bond.atom_b;
        else if (bond.atom_b == centre_atom) other = bond.atom_a;
        // A sizeable negative Mayer value is antibonding/non-structural
        // evidence, not a coordination contact.  Keep the sign contract in
        // sync with the molecular interaction graph and bond analyser.
        const double mayer = bond.mayer_order;
        if (other >= wavefunction.atoms.size() || other == centre_atom ||
            mayer < std::max(0.0, options.minimum_mayer_bond_order)) {
            continue;
        }
        const auto duplicate = std::find_if(
            bonded.begin(), bonded.end(),
            [&](const auto& item) { return item.first == other; });
        if (duplicate == bonded.end()) bonded.push_back({other, mayer});
        else duplicate->second = std::max(duplicate->second, mayer);
    }

    std::vector<CoordinationContact> candidates;
    candidates.reserve(bonded.size());
    for (const auto& [atom, mayer] : bonded) {
        auto contact = make_contact(wavefunction, centre_atom, atom, mayer);
        if (contact.distance > 1.0e-12) candidates.push_back(contact);
    }
    std::sort(candidates.begin(), candidates.end(), [](const auto& left,
                                                       const auto& right) {
        if (left.distance != right.distance) return left.distance < right.distance;
        return left.mayer_bond_order > right.mayer_bond_order;
    });

    std::vector<CoordinationContact> radial;
    for (const auto& candidate : candidates) {
        const bool shadowed = std::any_of(
            radial.begin(), radial.end(), [&](const auto& retained) {
                return direction_dot(candidate, retained) >
                       std::clamp(options.radial_shadow_cosine, -1.0, 1.0);
            });
        if (!shadowed) radial.push_back(candidate);
    }
    result.radial_candidate_count = radial.size();

    const double strongest = std::accumulate(
        radial.begin(), radial.end(), 0.0,
        [](const double value, const auto& item) {
            return std::max(value, item.mayer_bond_order);
        });
    const double electronic_floor = std::max(
        std::max(0.0, options.absolute_electronic_floor),
        std::max(0.0, options.relative_electronic_floor) * strongest);
    auto electronically_strong = radial;
    electronically_strong.erase(
        std::remove_if(electronically_strong.begin(),
                       electronically_strong.end(), [&](const auto& contact) {
                           return contact.mayer_bond_order < electronic_floor;
                       }),
        electronically_strong.end());

    CoordinationShell strong_shell;
    strong_shell.centre_atom = centre_atom;
    strong_shell.contacts = electronically_strong;
    strong_shell.radial_candidate_count = radial.size();
    const auto strong_match = analyse_coordination_shell(strong_shell, options.match);
    if (strong_match.accepted && !strong_match.ambiguous) {
        result.contacts = std::move(electronically_strong);
        return result;
    }

    CoordinationShell radial_shell;
    radial_shell.centre_atom = centre_atom;
    radial_shell.contacts = radial;
    radial_shell.radial_candidate_count = radial.size();
    const auto radial_match = analyse_coordination_shell(radial_shell, options.match);
    if (radial_match.accepted && !radial_match.ambiguous) {
        result.contacts = std::move(radial);
        result.used_low_mayer_retry = true;
        return result;
    }

    if (!electronically_strong.empty()) {
        result.contacts = std::move(electronically_strong);
    } else {
        result.contacts = std::move(radial);
        result.used_low_mayer_retry = true;
    }
    return result;
}

} // namespace cov
