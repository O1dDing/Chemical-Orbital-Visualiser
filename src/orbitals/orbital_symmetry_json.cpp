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

#include "cov/orbital_symmetry_scope.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include "cov/pi_pair_evidence.hpp"
#include <cmath>
#include <concepts>
#include <iomanip>
#include <sstream>

namespace cov {
namespace {
void quoted(std::ostream& out, const std::string& text) {
    out << '"';
    for (const unsigned char c : text) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u00" << "0123456789abcdef"[c >> 4] << "0123456789abcdef"[c & 15];
        else out << c;
    }
    out << '"';
}
void value(std::ostream& out, const std::string& v) { quoted(out, v); }
void value(std::ostream& out, const char* v) { quoted(out, v); }
void value(std::ostream& out, bool v) { out << (v ? "true" : "false"); }
template<std::integral T> void value(std::ostream& out, T v) { out << v; }
void value(std::ostream& out, double v) { if (std::isfinite(v)) out << v; else out << "null"; }
template<class T> void array(std::ostream& out, const T& values) {
    out << '['; bool first = true;
    for (const auto& v : values) { if (!first) out << ','; first = false; value(out, v); }
    out << ']';
}
struct Object {
    std::ostream& out;
    bool first = true;
    explicit Object(std::ostream& stream) : out(stream) { out << '{'; }
    ~Object() { out << '}'; }
    void key(const char* name) { if (!first) out << ','; first = false; quoted(out, name); out << ':'; }
    template<class T> void field(const char* name, const T& v) { key(name); value(out, v); }
    template<class T> void array_field(const char* name, const T& v) { key(name); array(out, v); }
};
const char* status(MetricSubspaceStatus state) {
    switch (state) {
        case MetricSubspaceStatus::Available: return "available";
        case MetricSubspaceStatus::InvalidInput: return "invalid-input";
        case MetricSubspaceStatus::Failed: return "failed";
        default: return "missing-input";
    }
}
void diagnostics(std::ostream& out, const MetricSubspaceDiagnostics& v) {
    Object o(out);
    o.field("supplied_columns", v.supplied_columns); o.field("numerical_rank", v.numerical_rank);
    o.field("unresolved_norm_columns", v.unresolved_norm_columns); o.field("rank_cutoff", v.rank_cutoff);
    o.field("gram_inverse_residual", v.gram_inverse_residual);
    o.field("orthonormality_residual", v.orthonormality_residual);
    o.field("retained_condition_number", v.retained_condition_number);
}
void metric(std::ostream& out, const MetricSubspaceOverlap& v) {
    Object o(out); o.field("status", status(v.status)); o.field("basis_dimension", v.basis_dimension);
    o.field("metric_numerical_rank", v.metric_numerical_rank);
    o.field("basis_norm_preconditioned", v.basis_norm_preconditioned);
    o.field("metric_rank_cutoff", v.metric_rank_cutoff); o.field("metric_symmetry_error", v.metric_symmetry_error);
    o.field("metric_eigen_residual", v.metric_eigen_residual);
    o.field("metric_minimum_eigenvalue", v.metric_minimum_eigenvalue);
    o.field("subspace_overlap_trace", v.subspace_overlap_trace); o.field("mean_orbital_fraction", v.mean_orbital_fraction);
    o.key("reference"); diagnostics(out, v.reference); o.key("orbitals"); diagnostics(out, v.orbitals);
    o.field("detail", v.detail);
}
void projection(std::ostream& out, const LocalAngularProjection& v) {
    Object o(out); o.field("status", status(v.status)); o.field("atom_index", v.atom_index);
    o.array_field("rotation_reference_to_input", v.rotation_reference_to_input);
    o.array_field("source_shell_indices", v.source_shell_indices);
    o.field("represented_spin_orbital_rank", v.represented_spin_orbital_rank);
    o.field("centre_projection_trace", v.centre_projection_trace); o.field("centre_mean_fraction", v.centre_mean_fraction);
    o.array_field("component_projection_traces", v.component_projection_traces);
    o.field("angular_partition_residual", v.angular_partition_residual); o.field("detail", v.detail);
    o.key("spins"); out << '['; bool first = true;
    for (const auto& s : v.spins) {
        if (!first) out << ','; first = false;
        Object spin(out); spin.field("spin", s.spin == Spin::Alpha ? "alpha" : "beta");
        spin.array_field("orbital_indices", s.orbital_indices); spin.key("centre"); metric(out, s.centre);
        spin.field("angular_partition_residual", s.angular_partition_residual);
        spin.key("components"); out << '['; bool first_component = true;
        for (const auto& c : s.components) {
            if (!first_component) out << ','; first_component = false;
            Object component(out); component.field("angular_degree", c.angular_degree);
            component.field("component", c.component); component.array_field("source_shell_indices", c.source_shell_indices);
            component.key("metric"); metric(out, c.metric);
            component.field("fraction_of_local_projection", c.fraction_of_local_projection);
        }
        out << ']';
    }
    out << ']';
}
void explanation(std::ostream& out, const OrbitalSymmetryExplanation& v) {
    Object o(out); o.field("label", v.label); o.field("origin", orbital_symmetry_origin_name(v.origin));
    o.field("point_group", v.point_group); o.field("point_group_basis", v.point_group_basis);
    o.array_field("orbital_indices", v.orbital_indices); o.array_field("atom_indices", v.atom_indices);
    o.field("axes_available", v.axes_available); o.key("rotation_reference_to_input");
    if (v.axes_available) array(out, v.rotation_reference_to_input); else out << "null";
    o.field("geometry_id", v.geometry_id); o.field("geometry_angular_rms", v.geometry_angular_rms);
    o.field("geometry_shape_measure", v.geometry_shape_measure); o.field("source_path", v.source_path);
    o.field("source_line_begin", v.source_line_begin); o.field("source_line_end", v.source_line_end);
    o.field("source_job_segment", v.source_job_segment);
    o.field("producer_detected_group", v.producer_detected_group);
    o.field("producer_abelian_group", v.producer_abelian_group); o.field("candidate_score", v.candidate_score);
    o.field("assessment_detail",v.assessment_detail);
    o.key("pi_partner_evidence");
    if(v.pi_partner_evidence)out<<pi_partner_assessment_json(*v.pi_partner_evidence);else out<<"null";
    o.key("local_decomposition");
    if(v.local_decomposition)projection(out,*v.local_decomposition);else out<<"null";
    o.key("candidate_source");
    if (v.candidate_source) explanation(out, *v.candidate_source); else out << "null";
    o.key("local_assignment");
    if (!v.local_assignment) out << "null";
    else {
        const auto& a = *v.local_assignment; Object item(out); item.field("label", a.label);
        item.field("source", a.source == LocalIrrepSource::MetricAngularProjection ? "metric-angular-projection" : "dimension-candidate");
        item.field("shell", static_cast<int>(a.shell)); item.field("copy_index", static_cast<int>(a.copy_index));
        item.field("basis_functions", a.basis_functions); item.field("conditional_shell_purity", a.confidence);
        item.field("centre_mean_fraction", a.centre_mean_fraction);
        item.field("angular_fraction_within_centre", a.angular_fraction_within_centre);
        item.field("labelled_fraction_of_target", a.labelled_fraction_of_target);
        item.field("measure", "unweighted-S-metric-subspace-overlap-not-electron-population");
        item.key("projection"); if (a.projection) projection(out, *a.projection); else out << "null";
    }
    o.key("molecular_assignment");
    if (!v.molecular_assignment) out << "null";
    else {
        const auto& a = *v.molecular_assignment; Object item(out);
        item.field("point_group", a.point_group); item.field("label", a.label);
        item.array_field("orbital_indices", a.orbital_indices); item.field("subspace_retention", a.subspace_retention);
        item.field("maximum_character_error", a.maximum_character_error);
        item.array_field("centre_bohr", a.centre_bohr); item.field("axes_available", a.axes_available);
        item.key("principal_axis"); if (a.axes_available) array(out, a.principal_axis); else out << "null";
        item.key("secondary_axis"); if (a.axes_available) array(out, a.secondary_axis); else out << "null";
        item.field("axis_convention", a.axis_convention);
    }
}
}
std::string orbital_symmetry_json(const OrbitalSymmetryExplanation& v) {
    std::ostringstream out; out << std::setprecision(std::numeric_limits<double>::max_digits10);
    explanation(out, v); return out.str();
}
}
