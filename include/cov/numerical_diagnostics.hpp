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

#include "cov/model.hpp"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <string_view>

namespace cov {
inline const char* numerical_status_name(const NumericalStatus status) noexcept {
    switch (status) {
        case NumericalStatus::Available:return "available";
        case NumericalStatus::MissingInput:return "missing-input";
        case NumericalStatus::InvalidInput:return "invalid-input";
        case NumericalStatus::Failed:return "failed";
        default:return "not-computed";
    }
}

namespace numerical_json {
inline void number(std::ostream& out,const double value) {
    if (std::isfinite(value)) out<<std::setprecision(17)<<value;
    else out<<"null";
}
inline void string(std::ostream& out,const std::string_view value) {
    out<<'"';
    constexpr char hex[]="0123456789abcdef";
    for (const unsigned char c:value) {
        if (c=='"' || c=='\\') out<<'\\'<<static_cast<char>(c);
        else if (c<32) out<<"\\u00"<<hex[c/16]<<hex[c%16];
        else out<<static_cast<char>(c);
    }
    out<<'"';
}
inline void metric(std::ostream& out,const AoMetricDiagnostics& data) {
    out<<"{\"status\":";string(out,numerical_status_name(data.status));
    out<<",\"dimension\":"<<data.dimension<<",\"finite\":"<<(data.finite?"true":"false")
       <<",\"positive_semidefinite\":"<<(data.positive_semidefinite?"true":"false")
       <<",\"symmetry_error\":";number(out,data.symmetry_error);
    out<<",\"minimum_eigenvalue\":";number(out,data.minimum_eigenvalue);
    out<<",\"maximum_eigenvalue\":";number(out,data.maximum_eigenvalue);
    out<<",\"condition_number\":";number(out,data.condition_number);
    out<<",\"eigen_residual\":";number(out,data.eigen_residual);
    out<<",\"numerical_rank\":"<<data.numerical_rank<<",\"numerical_rank_tolerance\":";
    number(out,data.numerical_rank_tolerance);
    out<<",\"detail\":";string(out,data.detail);
    out<<",\"spin_blocks\":[";
    bool first=true;
    for (const auto& block:data.orbital_blocks) {
        if (!first) out<<',';
        first=false;
        out<<"{\"spin\":";string(out,block.spin==Spin::Alpha?"alpha":"beta");
        out<<",\"orbital_count\":"<<block.orbital_count<<",\"status\":";
        string(out,numerical_status_name(block.status));
        out<<",\"maximum_diagonal_error\":";number(out,block.maximum_diagonal_error);
        out<<",\"maximum_off_diagonal_error\":";number(out,block.maximum_off_diagonal_error);
        out<<",\"orthonormality_tolerance\":";number(out,block.orthonormality_tolerance);
        out<<",\"detail\":";string(out,block.detail);out<<'}';
    }
    out<<"]}";
}
}

// Shared serialization of measured production data; it does not recompute or
// substitute independent science. Unknown/nonfinite numbers use JSON null.
inline void write_numerical_diagnostics_json(std::ostream& out,const Wavefunction& wf) {
    out<<"{\"schema\":1,\"scope\":\"full AO basis and every supplied spin MO block\","
         "\"units\":\"dimensionless inner products\",\"source\":\"analytic basis integrals\",\"basis\":";
    numerical_json::metric(out,wf.ao_metric_diagnostics);
    out<<",\"producer\":";numerical_json::metric(out,wf.producer_ao_metric_diagnostics);
    out<<",\"producer_basis_comparison\":{\"status\":";
    numerical_json::string(out,numerical_status_name(wf.producer_ao_overlap_basis_status));
    out<<",\"maximum_absolute_error\":";numerical_json::number(out,wf.producer_ao_overlap_basis_error);
    out<<",\"tolerance\":";numerical_json::number(out,wf.producer_ao_overlap_basis_tolerance);
    out<<"}}";
}

// Read-only serialization of the matrices and diagnostics used by production.
// The native validator writes one record per load; no alternate reconstruction
// or tolerance-dependent calculation is performed by this serializer.
inline void write_density_evidence_json(std::ostream& out,const Wavefunction& wf) {
    using namespace numerical_json;
    const auto matrix=[&](const std::vector<double>& packed,DataProvenance provenance,
                          const DensityMatrixDiagnostics& info) {
        out<<"{\"status\":";string(out,numerical_status_name(info.status));
        out<<",\"provenance\":";string(out,data_provenance_name(provenance));
        out<<",\"input_provenance\":";string(out,data_provenance_name(info.input_provenance));
        out<<",\"occupation_model\":";string(out,orbital_occupation_model_name(info.occupation_model));
        out<<",\"model_provenance\":";string(out,data_provenance_name(info.model_provenance));
        out<<",\"occupied_space_complete\":"<<(info.occupied_space_complete?"true":"false");
        out<<",\"occupation_sum\":";number(out,info.occupation_sum);
        out<<",\"expected_electron_count\":";number(out,info.expected_electron_count);
        out<<",\"metric_trace\":";number(out,info.metric_trace);
        out<<",\"reason\":";string(out,info.detail);
        out<<",\"packed\":[";
        for(std::size_t i=0;i<packed.size();++i){if(i)out<<',';number(out,packed[i]);}
        out<<"]}";
    };
    out<<"{\"schema\":1,\"basis_count\":"<<wf.basis_count
       <<",\"ao_representation\":\"COV internal normalized real or Cartesian AO order\","
         "\"packing\":\"lower triangular: i*(i+1)/2+j, i>=j\","
         "\"scalar_type\":\"IEEE754-float64\","
         "\"matrix_units\":\"electron occupation coefficients in the AO expansion\","
         "\"trace_units\":\"electrons\",\"source_format\":";
    string(out,wavefunction_source_name(wf.source));
    out<<",\"actual\":{\"total\":";
    matrix(wf.total_density_packed,wf.total_density_provenance,wf.total_density_diagnostics);
    out<<",\"spin\":";
    matrix(wf.spin_density_packed,wf.spin_density_provenance,wf.spin_density_diagnostics);
    out<<"},\"production_mo_reconstruction\":{\"total\":";
    matrix(wf.mo_density_reconstruction.total.packed,wf.mo_density_reconstruction.total.provenance,
           wf.mo_density_reconstruction.total.diagnostics);
    out<<",\"spin\":";
    matrix(wf.mo_density_reconstruction.spin.packed,wf.mo_density_reconstruction.spin.provenance,
           wf.mo_density_reconstruction.spin.diagnostics);
    out<<"}}";
}
}
