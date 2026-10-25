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

#include "cov/local_angular_projection.hpp"
#include "cov/local_angular_generators.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <utility>

namespace cov {
namespace {
std::vector<double> row_major_columns(std::size_t rows,const std::vector<std::vector<double>>& columns) {
    std::vector<double> matrix(rows*columns.size());
    for(std::size_t j=0;j<columns.size();++j)for(std::size_t i=0;i<rows;++i)
        matrix[i*columns.size()+j]=columns[j][i];
    return matrix;
}
struct ReferenceComponent {
    int degree=0;
    std::size_t component=0;
    std::vector<std::size_t> shells;
    MetricPreparedSubspace space;
};
}
struct LocalAngularProjectionWorkspace::Impl {
    MetricSubspaceStatus status=MetricSubspaceStatus::MissingInput;
    std::string detail;
    std::size_t dimension=0,atom=0;
    std::array<double,9> rotation{};
    std::vector<std::size_t> shells;
    std::vector<std::vector<double>> orbitals;
    std::vector<Spin> spins;
    std::unique_ptr<MetricSubspaceContext> metric;
    MetricPreparedSubspace centre;
    std::array<ReferenceComponent,25> components;
};

LocalAngularProjectionWorkspace::LocalAngularProjectionWorkspace(const Wavefunction& wf,std::size_t atom,
    const std::array<double,9>& rotation) {
    auto data=std::make_shared<Impl>();
    impl_=data;
    data->atom=atom;data->rotation=rotation;data->dimension=wf.basis_count;
    auto invalid=[&](const char* detail){data->status=MetricSubspaceStatus::InvalidInput;data->detail=detail;};
    if(atom>=wf.atoms.size()) {invalid("Requested local centre is outside the input atoms");return;}
    if(wf.basis_count==0 || wf.ao_overlap.empty() || wf.shells.empty()) {
        data->detail="Local angular projection requires AO definitions and the AO overlap metric";return;
    }
    const std::size_t n=wf.basis_count;
    if(n>std::numeric_limits<std::size_t>::max()/n || wf.ao_overlap.size()!=n*n) {
        invalid("AO metric dimensions do not match the input basis");return;
    }
    data->metric=std::make_unique<MetricSubspaceContext>(n,wf.ao_overlap);
    const auto metric_info=data->metric->diagnostics();
    if(metric_info.status!=MetricSubspaceStatus::Available) {
        data->status=metric_info.status;data->detail=metric_info.detail;return;
    }
    std::vector<bool> covered(n,false);
    std::vector<std::vector<double>> centre_columns;
    std::array<std::vector<std::vector<double>>,25> component_columns;
    for(int l=0;l<=4;++l)for(std::size_t m=0;m<static_cast<std::size_t>(2*l+1);++m) {
        auto& ref=data->components[static_cast<std::size_t>(l*l)+m];ref.degree=l;ref.component=m;
    }
    for(std::size_t shell_index=0;shell_index<wf.shells.size();++shell_index) {
        const auto& shell=wf.shells[shell_index];
        const std::size_t count=static_cast<std::size_t>(shell_basis_count(shell));
        const std::size_t offset=shell.basis_offset;
        if(shell.atom_index>=wf.atoms.size() || shell.angular_momentum>4 || offset>n || count>n-offset) {
            invalid("AO shell range or supported angular degree is invalid");return;
        }
        for(std::size_t mu=offset;mu<offset+count;++mu) {
            if(covered[mu]) {invalid("AO shell ranges overlap");return;}
            covered[mu]=true;
        }
        if(shell.atom_index!=atom)continue;
        data->shells.push_back(shell_index);
        for(std::size_t mu=offset;mu<offset+count;++mu) {
            std::vector<double> v(n,0);v[mu]=1;centre_columns.push_back(std::move(v));
        }
        for(int l=shell.angular_momentum;l>=0;l-=2) {
            if(shell.pure && l!=shell.angular_momentum)continue;
            const auto generators=local_angular_generators(shell.angular_momentum,shell.pure!=0,l,rotation);
            if(!generators.available || generators.basis_count!=count) {
                invalid("Local angular reference cannot be represented in the supplied AO frame");return;
            }
            for(std::size_t component=0;component<generators.columns;++component) {
                const auto index=static_cast<std::size_t>(l*l)+component;
                std::vector<double> column(n,0);
                for(std::size_t mu=0;mu<count;++mu)
                    column[offset+mu]=generators.coefficients[mu*generators.columns+component];
                component_columns[index].push_back(std::move(column));
                data->components[index].shells.push_back(shell_index);
            }
        }
    }
    if(std::find(covered.begin(),covered.end(),false)!=covered.end()) {
        invalid("Declared AO basis has an unrepresented direction");return;
    }
    if(centre_columns.empty()) {data->detail="The requested atom has no represented AO shell";return;}
    data->centre=data->metric->prepare_space(centre_columns.size(),row_major_columns(n,centre_columns));
    if(data->centre.status()!=MetricSubspaceStatus::Available) {
        data->status=data->centre.status();data->detail="Local centre reference preparation failed";return;
    }
    for(std::size_t i=0;i<25;++i) {
        data->components[i].space=data->metric->prepare_space(component_columns[i].size(),row_major_columns(n,component_columns[i]));
        if(data->components[i].space.status()!=MetricSubspaceStatus::Available) {
            data->status=data->components[i].space.status();data->detail="Local angular reference preparation failed";return;
        }
    }
    for(const auto& mo:wf.orbitals) {
        if(mo.coefficients.size()!=n || (mo.spin!=Spin::Alpha && mo.spin!=Spin::Beta)) {
            invalid("Orbital dimensions or spin identity are invalid");return;
        }
        data->orbitals.push_back(mo.coefficients);data->spins.push_back(mo.spin);
    }
    data->status=MetricSubspaceStatus::Available;
    data->detail="Unweighted S-metric local angular subspace projections; not additive atom populations or a whole-molecule irrep assignment";
}

LocalAngularProjection LocalAngularProjectionWorkspace::project(std::span<const std::size_t> indices) const {
    LocalAngularProjection result;
    const auto& data=*impl_;
    result.status=data.status;result.atom_index=data.atom;result.rotation_reference_to_input=data.rotation;
    result.source_shell_indices=data.shells;result.detail=data.detail;
    if(data.status!=MetricSubspaceStatus::Available)return result;
    if(indices.empty()) {result.status=MetricSubspaceStatus::MissingInput;result.detail="No target orbital subspace supplied";return result;}
    std::set<std::size_t> seen;
    std::array<std::vector<std::size_t>,2> by_spin;
    for(auto index:indices) {
        if(index>=data.orbitals.size() || !seen.insert(index).second) {
            result.status=MetricSubspaceStatus::InvalidInput;result.detail="Invalid or repeated target orbital identity";return result;
        }
        by_spin[data.spins[index]==Spin::Beta?1:0].push_back(index);
    }
    auto failed=[&](MetricSubspaceStatus status,const std::string& detail) {
        result.status=status;result.detail=detail;
        result.centre_projection_trace=std::numeric_limits<double>::quiet_NaN();
        result.centre_mean_fraction=std::numeric_limits<double>::quiet_NaN();
        result.component_projection_traces.fill(std::numeric_limits<double>::quiet_NaN());
    };
    result.centre_projection_trace=0;
    result.component_projection_traces.fill(0);
    for(std::size_t spin=0;spin<2;++spin) {
        if(by_spin[spin].empty())continue;
        LocalAngularSpinProjection block;
        block.spin=spin==1?Spin::Beta:Spin::Alpha;block.orbital_indices=by_spin[spin];
        std::vector<std::vector<double>> columns;
        for(auto index:block.orbital_indices)columns.push_back(data.orbitals[index]);
        const auto target=data.metric->prepare_space(columns.size(),row_major_columns(data.dimension,columns));
        block.centre=data.metric->compare(data.centre,target);
        if(block.centre.status!=MetricSubspaceStatus::Available) {
            failed(block.centre.status,block.centre.detail);return result;
        }
        double component_sum=0;
        for(std::size_t i=0;i<25;++i) {
            const auto& reference=data.components[i];
            LocalAngularProjectionComponent component;
            component.angular_degree=reference.degree;component.component=reference.component;
            component.source_shell_indices=reference.shells;
            component.metric=data.metric->compare(reference.space,target);
            if(component.metric.status!=MetricSubspaceStatus::Available) {
                failed(component.metric.status,component.metric.detail);return result;
            }
            const double trace=component.metric.subspace_overlap_trace;
            component_sum+=trace;result.component_projection_traces[i]+=trace;
            if(block.centre.subspace_overlap_trace>0)
                component.fraction_of_local_projection=trace/block.centre.subspace_overlap_trace;
            block.components.push_back(std::move(component));
        }
        block.angular_partition_residual=component_sum-block.centre.subspace_overlap_trace;
        result.centre_projection_trace+=block.centre.subspace_overlap_trace;
        result.represented_spin_orbital_rank+=block.centre.orbitals.numerical_rank;
        result.spins.push_back(std::move(block));
    }
    double component_sum=0;for(double trace:result.component_projection_traces)component_sum+=trace;
    result.angular_partition_residual=component_sum-result.centre_projection_trace;
    if(result.represented_spin_orbital_rank>0)
        result.centre_mean_fraction=result.centre_projection_trace/static_cast<double>(result.represented_spin_orbital_rank);
    return result;
}
}
