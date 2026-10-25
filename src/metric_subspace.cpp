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

#include "cov/metric_subspace.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cov {
namespace {
using Matrix=Eigen::MatrixXd;
using RowMatrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
constexpr double eps=std::numeric_limits<double>::epsilon();
bool valid_size(std::size_t rows,std::size_t columns,std::size_t size) {
    return rows && columns<=std::numeric_limits<std::size_t>::max()/rows && rows*columns==size;
}
struct Prepared { Matrix columns; MetricSubspaceDiagnostics diagnostics; };
Matrix balance_coefficients(const Matrix& input,const Eigen::VectorXd& root_norms) {
    Matrix result(input.rows(),input.cols());
    std::vector<double> mantissas(static_cast<std::size_t>(input.rows()));
    std::vector<int> exponents(static_cast<std::size_t>(input.rows()));
    for(Eigen::Index j=0;j<input.cols();++j) {
        int largest=std::numeric_limits<int>::min();
        for(Eigen::Index i=0;i<input.rows();++i) {
            int coefficient_exponent=0,norm_exponent=0;
            const double coefficient=std::frexp(input(i,j),&coefficient_exponent);
            const double norm=std::frexp(root_norms[i],&norm_exponent);
            const auto row=static_cast<std::size_t>(i);
            mantissas[row]=coefficient*norm;exponents[row]=coefficient_exponent+norm_exponent;
            if(coefficient!=0)largest=std::max(largest,exponents[row]);
        }
        for(Eigen::Index i=0;i<input.rows();++i) {
            const auto row=static_cast<std::size_t>(i);
            result(i,j)=mantissas[row]==0?0:std::ldexp(mantissas[row],exponents[row]-largest);
        }
    }
    return result;
}
Prepared prepare(const Matrix& metric,const Matrix& input) {
    Prepared result;auto& diagnostics=result.diagnostics;
    diagnostics.supplied_columns=static_cast<std::size_t>(input.cols());
    result.columns.resize(input.rows(),0);
    if(input.cols()==0)return result;
    Matrix scaled=input;
    for(Eigen::Index j=0;j<scaled.cols();++j) {
        const double maximum=scaled.col(j).cwiseAbs().maxCoeff();
        if(maximum>0)scaled.col(j)/=maximum;
    }
    const Matrix product=metric*scaled;
    const Matrix gram_raw=scaled.transpose()*product;
    Matrix normalized=scaled;
    for(Eigen::Index j=0;j<scaled.cols();++j) {
        const double norm=gram_raw(j,j);
        const double error=64*static_cast<double>(metric.rows())*eps*
            metric.cwiseAbs().maxCoeff()*scaled.col(j).cwiseAbs().sum()*scaled.col(j).cwiseAbs().sum();
        if(!std::isfinite(norm) || norm < -error)throw std::invalid_argument("Invalid subspace norm in AO metric");
        if(norm<=error) {
            if(scaled.col(j).squaredNorm()>0)++diagnostics.unresolved_norm_columns;
            normalized.col(j).setZero();
        } else normalized.col(j)/=std::sqrt(norm);
    }
    Matrix gram=normalized.transpose()*metric*normalized;
    gram=(gram+gram.transpose()).eval()*.5;
    Eigen::SelfAdjointEigenSolver<Matrix> solver(gram);
    if(solver.info()!=Eigen::Success)throw std::runtime_error("Subspace Gram eigensolver did not converge");
    const auto& values=solver.eigenvalues();
    const double scale=std::max(values.cwiseAbs().maxCoeff(),std::numeric_limits<double>::min());
    diagnostics.rank_cutoff=64*static_cast<double>(std::max(input.rows(),input.cols()))*eps*scale;
    if(values[0]<-diagnostics.rank_cutoff)throw std::invalid_argument("Indefinite subspace Gram matrix");
    diagnostics.numerical_rank=static_cast<std::size_t>((values.array()>diagnostics.rank_cutoff).count());
    const Eigen::Index rank=static_cast<Eigen::Index>(diagnostics.numerical_rank);
    Matrix inverse=Matrix::Zero(gram.rows(),gram.cols());
    result.columns.resize(input.rows(),rank);
    Eigen::Index column=0;
    for(Eigen::Index i=0;i<values.size();++i)if(values[i]>diagnostics.rank_cutoff) {
        const auto direction=solver.eigenvectors().col(i);
        result.columns.col(column++)=normalized*direction/std::sqrt(values[i]);
        inverse.noalias()+=(direction*direction.transpose())/values[i];
    }
    diagnostics.gram_inverse_residual=(gram*inverse*gram-gram).cwiseAbs().maxCoeff()/scale;
    if(rank) {
        diagnostics.orthonormality_residual=(result.columns.transpose()*metric*result.columns-Matrix::Identity(rank,rank)).cwiseAbs().maxCoeff();
        diagnostics.retained_condition_number=values[values.size()-1]/values[values.size()-rank];
    }
    if(!result.columns.allFinite() || !std::isfinite(diagnostics.gram_inverse_residual) || !std::isfinite(diagnostics.orthonormality_residual))
        throw std::runtime_error("Nonfinite metric subspace preparation");
    return result;
}
}


struct MetricSubspaceContext::Impl {
    Matrix metric;
    Eigen::VectorXd root_norms;
    MetricSubspaceOverlap report;
};
struct MetricPreparedSubspace::Impl {
    Prepared space;
    std::shared_ptr<const void> context_token;
    MetricSubspaceStatus status=MetricSubspaceStatus::MissingInput;
    std::string detail;
};

MetricSubspaceContext::MetricSubspaceContext(std::size_t dimension,const std::vector<double>& values) {
    auto data=std::make_shared<Impl>();impl_=data;
    auto& result=data->report;result.basis_dimension=dimension;
    if(!dimension || values.empty()) {result.detail="AO metric input is absent";return;}
    if(!valid_size(dimension,dimension,values.size())) {
        result.status=MetricSubspaceStatus::InvalidInput;result.detail="Incompatible AO metric dimensions";return;
    }
    try {
        const auto n=static_cast<Eigen::Index>(dimension);
        const double size=static_cast<double>(dimension);
        const Eigen::Map<const RowMatrix> supplied(values.data(),n,n);
        if(!supplied.allFinite())throw std::invalid_argument("Nonfinite AO metric input");
        data->root_norms.resize(n);
        for(Eigen::Index i=0;i<n;++i) {
            if(supplied(i,i)<=0)throw std::invalid_argument("AO basis function has no positive norm");
            data->root_norms[i]=std::sqrt(supplied(i,i));
        }
        Matrix metric(n,n);
        for(Eigen::Index i=0;i<n;++i)for(Eigen::Index j=0;j<n;++j)
            metric(i,j)=i==j?1.0:(supplied(i,j)/data->root_norms[i])/data->root_norms[j];
        if(!metric.allFinite())throw std::invalid_argument("Nonfinite norm-balanced AO metric");
        result.basis_norm_preconditioned=true;
        result.metric_symmetry_error=(metric-metric.transpose()).cwiseAbs().maxCoeff();
        if(result.metric_symmetry_error>64*size*eps)throw std::invalid_argument("Nonsymmetric AO metric");
        data->metric=(metric+metric.transpose())*.5;
        Eigen::SelfAdjointEigenSolver<Matrix> solver(data->metric);
        if(solver.info()!=Eigen::Success)throw std::runtime_error("AO metric eigensolver did not converge");
        const auto& eigenvalues=solver.eigenvalues();
        result.metric_rank_cutoff=64*size*eps*std::max(eigenvalues.cwiseAbs().maxCoeff(),std::numeric_limits<double>::min());
        result.metric_minimum_eigenvalue=eigenvalues[0];
        result.metric_numerical_rank=static_cast<std::size_t>((eigenvalues.array()>result.metric_rank_cutoff).count());
        result.metric_eigen_residual=(data->metric*solver.eigenvectors()-solver.eigenvectors()*eigenvalues.asDiagonal()).cwiseAbs().maxCoeff();
        if(eigenvalues[0]<-result.metric_rank_cutoff)throw std::invalid_argument("Indefinite AO metric");
        if(result.metric_eigen_residual>64*size*eps)throw std::runtime_error("AO eigensolver residual exceeded the arithmetic bound");
        result.status=MetricSubspaceStatus::Available;
        result.detail="Immutable AO-norm-balanced metric context; full spectrum inspected without truncation";
    } catch(const std::invalid_argument& error) {result.status=MetricSubspaceStatus::InvalidInput;result.detail=error.what();}
      catch(const std::exception& error) {result.status=MetricSubspaceStatus::Failed;result.detail=error.what();}
}

MetricSubspaceOverlap MetricSubspaceContext::diagnostics() const {return impl_->report;}

MetricPreparedSubspace MetricSubspaceContext::prepare_space(std::size_t columns,const std::vector<double>& coefficients) const {
    MetricPreparedSubspace result;auto data=std::make_shared<MetricPreparedSubspace::Impl>();result.impl_=data;
    data->context_token=impl_;data->space.diagnostics.supplied_columns=columns;
    data->status=impl_->report.status;data->detail=impl_->report.detail;
    if(data->status!=MetricSubspaceStatus::Available)return result;
    if(!valid_size(impl_->report.basis_dimension,columns,coefficients.size())) {
        data->status=MetricSubspaceStatus::InvalidInput;data->detail="Incompatible subspace matrix dimensions";return result;
    }
    try {
        const Eigen::Map<const RowMatrix> input(coefficients.data(),impl_->metric.rows(),static_cast<Eigen::Index>(columns));
        if(!input.allFinite())throw std::invalid_argument("Nonfinite subspace coefficients");
        data->space=prepare(impl_->metric,balance_coefficients(input,impl_->root_norms));
        data->status=MetricSubspaceStatus::Available;data->detail="Metric subspace prepared with explicit numerical rank and residuals";
    } catch(const std::invalid_argument& error) {data->status=MetricSubspaceStatus::InvalidInput;data->detail=error.what();}
      catch(const std::exception& error) {data->status=MetricSubspaceStatus::Failed;data->detail=error.what();}
    return result;
}

MetricSubspaceDiagnostics MetricPreparedSubspace::diagnostics() const {
    return impl_?impl_->space.diagnostics:MetricSubspaceDiagnostics{};
}
MetricSubspaceStatus MetricPreparedSubspace::status() const {
    return impl_?impl_->status:MetricSubspaceStatus::MissingInput;
}

MetricSubspaceOverlap MetricSubspaceContext::compare(
    const MetricPreparedSubspace& reference,const MetricPreparedSubspace& orbitals) const {
    auto result=impl_->report;
    if(result.status!=MetricSubspaceStatus::Available)return result;
    if(!reference.impl_ || !orbitals.impl_) {
        result.status=MetricSubspaceStatus::MissingInput;result.detail="Prepared subspace input is absent";return result;
    }
    if(reference.impl_->context_token.get()!=impl_.get() || orbitals.impl_->context_token.get()!=impl_.get()) {
        result.status=MetricSubspaceStatus::InvalidInput;result.detail="Prepared spaces belong to a different AO metric context";return result;
    }
    result.reference=reference.impl_->space.diagnostics;result.orbitals=orbitals.impl_->space.diagnostics;
    for(const auto* space:{reference.impl_.get(),orbitals.impl_.get()}) {
        if(space->status!=MetricSubspaceStatus::Available) {result.status=space->status;result.detail=space->detail;return result;}
    }
    result.subspace_overlap_trace=(reference.impl_->space.columns.transpose()*impl_->metric*orbitals.impl_->space.columns).squaredNorm();
    if(!std::isfinite(result.subspace_overlap_trace)) {
        result.status=MetricSubspaceStatus::Failed;result.detail="Nonfinite metric subspace overlap";return result;
    }
    if(result.orbitals.numerical_rank)result.mean_orbital_fraction=result.subspace_overlap_trace/static_cast<double>(result.orbitals.numerical_rank);
    result.detail=result.orbitals.numerical_rank?
        "S-metric projector overlap in one immutable AO-norm-balanced context; dimensionless trace; ranks and residuals retained":
        "Target has zero numerically resolved metric rank; subspace trace is zero and normalized mean is undefined";
    return result;
}

MetricSubspaceOverlap metric_subspace_overlap(std::size_t dimension,const std::vector<double>& metric,
    std::size_t reference_columns,const std::vector<double>& reference,
    std::size_t orbital_columns,const std::vector<double>& orbitals) {
    const MetricSubspaceContext context(dimension,metric);
    return context.compare(context.prepare_space(reference_columns,reference),context.prepare_space(orbital_columns,orbitals));
}
}
