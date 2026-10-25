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

#include "cov/density.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cov {
namespace {
struct OccupationWeights {
    std::vector<double> total, spin;
    DensityMatrixDiagnostics total_info, spin_info;
};
void fail(DensityMatrixDiagnostics& info, NumericalStatus status, const char* reason) {
    info.status=status;
    info.occupied_space_complete=false;
    info.detail=reason;
}
OccupationWeights occupation_weights(const Wavefunction& wf) {
    OccupationWeights result;
    auto& total=result.total_info; auto& spin=result.spin_info;
    total.input_provenance=spin.input_provenance=DataProvenance::Derived;
    const auto fail_both=[&](NumericalStatus status,const char* reason) {
        fail(total,status,reason); fail(spin,status,reason);
    };
    if(wf.basis_count==0) {
        fail_both(NumericalStatus::MissingInput,"AO basis dimension is missing"); return result;
    }
    const bool fchk=wf.source==WavefunctionSource::Fchk;
    const bool canonical=wf.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared;
    const bool explicit_beta=std::any_of(wf.orbitals.begin(),wf.orbitals.end(),[](const auto& mo) {
        return mo.spin_provenance!=DataProvenance::Unavailable && mo.spin==Spin::Beta;
    });
    const bool explicit_spin=wf.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin ||
                             (!canonical && explicit_beta);
    const bool shared_model=canonical || wf.source==WavefunctionSource::Molden ||
        wf.orbital_occupation_model==OrbitalOccupationModel::SharedIntegerDeterminant ||
        wf.orbital_occupation_model==OrbitalOccupationModel::SharedFractionalUnresolved;
    // An inferred shared determinant is named only after every occupation and
    // spin field has been checked. A default enum is not source evidence.
    const auto model=canonical?OrbitalOccupationModel::CanonicalShared:
        explicit_spin?OrbitalOccupationModel::ExplicitSpin:OrbitalOccupationModel::Unspecified;
    total.occupation_model=spin.occupation_model=model;
    total.model_provenance=spin.model_provenance=model==OrbitalOccupationModel::Unspecified?
        DataProvenance::Unavailable:
        wf.orbital_occupation_model==model && wf.orbital_occupation_model_provenance!=DataProvenance::Unavailable?
            wf.orbital_occupation_model_provenance:DataProvenance::Derived;
    if(fchk && wf.electron_counts_provenance==DataProvenance::Unavailable) {
        fail_both(NumericalStatus::MissingInput,"Canonical FCHK electron counts are unavailable"); return result;
    }
    if(fchk) {
        total.expected_electron_count=static_cast<double>(wf.alpha_electrons)+wf.beta_electrons;
        spin.expected_electron_count=static_cast<double>(wf.alpha_electrons)-wf.beta_electrons;
        if(wf.alpha_electrons>wf.basis_count || wf.beta_electrons>wf.basis_count) {
            fail_both(NumericalStatus::InvalidInput,"Canonical occupied count exceeds the AO dimension"); return result;
        }
    }
    if(wf.orbitals.empty() && !(fchk && total.expected_electron_count==0.0)) {
        fail_both(NumericalStatus::MissingInput,"MO occupations and coefficient directions are missing"); return result;
    }
    bool spins_known=true, integer_occupations=true;
    double total_count=0.0,spin_count=0.0,alpha_count=0.0,beta_count=0.0;
    std::unordered_set<std::size_t> alpha_indices,beta_indices;
    result.total.reserve(wf.orbitals.size());result.spin.reserve(wf.orbitals.size());
    for(const auto& mo:wf.orbitals) {
        if(mo.occupation_provenance==DataProvenance::Unavailable) {
            fail_both(NumericalStatus::MissingInput,"At least one MO Occup field is missing"); return result;
        }
        const double occupation=mo.occupation;
        const double maximum=explicit_spin?1.0:2.0;
        if(!std::isfinite(occupation) || occupation<0.0 || occupation>maximum) {
            fail_both(NumericalStatus::InvalidInput,"MO occupation is nonfinite or outside its model range"); return result;
        }
        const bool known_spin=mo.spin_provenance!=DataProvenance::Unavailable;
        spins_known=spins_known && known_spin;
        integer_occupations=integer_occupations && std::floor(occupation)==occupation;
        result.total.push_back(occupation);total_count+=occupation;
        double spin_weight=0.0;
        if(fchk && (canonical || explicit_spin)) {
            if(!known_spin || mo.source_orbital_index==std::numeric_limits<std::size_t>::max()) {
                fail_both(NumericalStatus::MissingInput,"Canonical source MO identity is missing"); return result;
            }
            if(mo.source_orbital_index>=wf.basis_count || (canonical && mo.spin!=Spin::Alpha)) {
                fail_both(NumericalStatus::InvalidInput,"Canonical source MO identity is inconsistent"); return result;
            }
            auto& indices=mo.spin==Spin::Beta?beta_indices:alpha_indices;
            if(!indices.insert(mo.source_orbital_index).second) {
                fail_both(NumericalStatus::InvalidInput,"A canonical source MO identity is duplicated"); return result;
            }
            const double alpha=mo.source_orbital_index<wf.alpha_electrons?1.0:0.0;
            const double beta=mo.source_orbital_index<wf.beta_electrons?1.0:0.0;
            const double expected=canonical?alpha+beta:mo.spin==Spin::Alpha?alpha:beta;
            if(occupation!=expected) {
                fail_both(NumericalStatus::InvalidInput,"MO occupation contradicts its canonical source identity"); return result;
            }
            spin_weight=canonical?alpha-beta:mo.spin==Spin::Alpha?occupation:-occupation;
        } else if(explicit_spin) {
            spin_weight=mo.spin==Spin::Beta?-occupation:occupation;
        } else if(shared_model && std::floor(occupation)==occupation) {
            // Explicit conventional shared determinant assumption; no index refill.
            spin_weight=occupation==1.0?1.0:0.0;
        }
        result.spin.push_back(spin_weight);spin_count+=spin_weight;
        alpha_count+=(occupation+spin_weight)*0.5;beta_count+=(occupation-spin_weight)*0.5;
    }
    if(!std::isfinite(total_count) || !std::isfinite(spin_count)) {
        fail_both(NumericalStatus::Failed,"Occupation accumulation overflowed");return result;
    }
    total.occupation_sum=total_count;spin.occupation_sum=spin_count;
    if(!canonical && !explicit_spin && shared_model && spins_known) {
        total.occupation_model=spin.occupation_model=integer_occupations?
            OrbitalOccupationModel::SharedIntegerDeterminant:OrbitalOccupationModel::SharedFractionalUnresolved;
        total.model_provenance=spin.model_provenance=DataProvenance::Derived;
    }
    total.status=spin.status=NumericalStatus::Available;
    total.occupied_space_complete=spin.occupied_space_complete=true;
    total.detail="Sum of all supplied MO occupations; completeness is relative to the declared input model";
    spin.detail=explicit_spin?"Signed occupations of explicit alpha and beta orbitals":
        canonical?"Alpha-minus-beta weights follow immutable canonical source MO identities":
        "Conventional shared integer determinant: Occup=2 paired, Occup=1 alpha, Occup=0 empty";
    if(fchk && !canonical && !explicit_spin) {
        fail_both(NumericalStatus::MissingInput,"FCHK source orbital layout was not established");
    } else if(fchk) {
        if(total_count!=total.expected_electron_count)
            fail(total,NumericalStatus::MissingInput,"Canonical occupied MO space is truncated");
        // Equal counts do not make different unrestricted occupied spaces cancel.
        // In a shared closed shell, all spin weights are explicitly zero.
        if((explicit_spin && (alpha_count!=wf.alpha_electrons || beta_count!=wf.beta_electrons)) ||
           (canonical && spin_count!=spin.expected_electron_count))
            fail(spin,NumericalStatus::MissingInput,"Spin-bearing canonical occupied MO directions are missing");
    } else if(!spins_known) {
        spin.occupation_sum=std::numeric_limits<double>::quiet_NaN();
        fail(spin,NumericalStatus::MissingInput,"At least one MO Spin field is missing or unrecognized");
    } else if(!explicit_spin && !shared_model) {
        spin.occupation_sum=std::numeric_limits<double>::quiet_NaN();
        fail(spin,NumericalStatus::MissingInput,"No spin partition or shared spatial occupation model was declared");
    } else if(!explicit_spin && !integer_occupations) {
        spin.occupation_sum=std::numeric_limits<double>::quiet_NaN();
        fail(spin,NumericalStatus::MissingInput,"Shared fractional occupations do not determine an alpha/beta partition");
    }
    return result;
}
DensityMatrixResult matrix_from_weights(const Wavefunction& wf,const std::vector<double>& weights,
                                       DensityMatrixDiagnostics info) {
    DensityMatrixResult result;result.diagnostics=std::move(info);
    if(result.diagnostics.status!=NumericalStatus::Available)return result;
    if(weights.size()!=wf.orbitals.size()) {
        fail(result.diagnostics,NumericalStatus::Failed,"Internal occupation-weight dimension mismatch");return result;
    }
    const std::size_t n=wf.basis_count;
    result.packed.assign(n*(n+1u)/2u,0.0);
    for(std::size_t orbital=0;orbital<weights.size();++orbital) {
        const double weight=weights[orbital];if(weight==0.0)continue;
        const auto& coefficients=wf.orbitals[orbital].coefficients;
        if(coefficients.size()!=n) {
            result.packed.clear();
            fail(result.diagnostics,NumericalStatus::MissingInput,"A required occupied MO coefficient direction is incomplete");return result;
        }
        if(!std::all_of(coefficients.begin(),coefficients.end(),[](double x){return std::isfinite(x);})) {
            result.packed.clear();
            fail(result.diagnostics,NumericalStatus::InvalidInput,"A required occupied MO coefficient is nonfinite");return result;
        }
        for(std::size_t i=0;i<n;++i) {
            const auto row=i*(i+1u)/2u;
            for(std::size_t j=0;j<=i;++j)
                result.packed[row+j]+=weight*coefficients[i]*coefficients[j];
        }
    }
    if(!std::all_of(result.packed.begin(),result.packed.end(),[](double x){return std::isfinite(x);})) {
        result.packed.clear();
        fail(result.diagnostics,NumericalStatus::Failed,"AO density accumulation produced a nonfinite value");return result;
    }
    result.provenance=DataProvenance::Derived;
    return result;
}
void establish_component(const Wavefunction& wf,std::vector<double>& actual,DataProvenance& provenance,
                         DensityMatrixDiagnostics& info,const DensityMatrixResult& reconstructed,
                         bool spin,bool reconstruct_missing) {
    if(provenance==DataProvenance::Producer) {
        info={};info.input_provenance=DataProvenance::Producer;
        const std::size_t n=wf.basis_count;
        if(n==0 || actual.size()!=n*(n+1u)/2u ||
           !std::all_of(actual.begin(),actual.end(),[](double x){return std::isfinite(x);})) {
            actual.clear();provenance=DataProvenance::Unavailable;
            fail(info,NumericalStatus::InvalidInput,"Producer AO density has an invalid dimension or nonfinite entry");return;
        }
        info.status=NumericalStatus::Available;info.occupied_space_complete=true;
        info.detail="Complete producer AO matrix; independent of supplied occupied MO coverage";
        if(wf.electron_counts_provenance!=DataProvenance::Unavailable)
            info.expected_electron_count=static_cast<double>(wf.alpha_electrons)+
                (spin?-static_cast<double>(wf.beta_electrons):static_cast<double>(wf.beta_electrons));
    } else if(info.input_provenance==DataProvenance::Producer && info.status==NumericalStatus::InvalidInput) {
        // Repeated processing cannot conceal a rejected producer matrix by fallback.
        actual.clear();provenance=DataProvenance::Unavailable;
    } else if(reconstruct_missing) {
        actual=reconstructed.packed;provenance=reconstructed.provenance;info=reconstructed.diagnostics;
    } else {
        actual.clear();provenance=DataProvenance::Unavailable;info={};
        fail(info,NumericalStatus::MissingInput,"No producer density; MO reconstruction was disabled");
    }
}
double metric_trace(const Wavefunction& wf,const std::vector<double>& packed) {
    const std::size_t n=wf.basis_count;
    if(n==0 || packed.size()!=n*(n+1u)/2u || wf.ao_overlap.size()!=n*n)
        return std::numeric_limits<double>::quiet_NaN();
    double trace=0.0;
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<=i;++j) {
        const double p=packed[i*(i+1u)/2u+j];
        trace+=p*wf.ao_overlap[j*n+i];
        if(i!=j)trace+=p*wf.ao_overlap[i*n+j];
    }
    return std::isfinite(trace)?trace:std::numeric_limits<double>::quiet_NaN();
}
} // namespace

DensityReconstruction reconstruct_density(const Wavefunction& wf) {
    auto weights=occupation_weights(wf);
    return {matrix_from_weights(wf,weights.total,std::move(weights.total_info)),
            matrix_from_weights(wf,weights.spin,std::move(weights.spin_info))};
}
void establish_density(Wavefunction& wf,bool reconstruct_missing) {
    wf.mo_density_reconstruction=reconstruct_density(wf);
    establish_component(wf,wf.total_density_packed,wf.total_density_provenance,wf.total_density_diagnostics,
                        wf.mo_density_reconstruction.total,false,reconstruct_missing);
    establish_component(wf,wf.spin_density_packed,wf.spin_density_provenance,wf.spin_density_diagnostics,
                        wf.mo_density_reconstruction.spin,true,reconstruct_missing);
    update_density_metric_diagnostics(wf);
}
void update_density_metric_diagnostics(Wavefunction& wf) {
    wf.total_density_diagnostics.metric_trace=metric_trace(wf,wf.total_density_packed);
    wf.spin_density_diagnostics.metric_trace=metric_trace(wf,wf.spin_density_packed);
    wf.mo_density_reconstruction.total.diagnostics.metric_trace=metric_trace(wf,wf.mo_density_reconstruction.total.packed);
    wf.mo_density_reconstruction.spin.diagnostics.metric_trace=metric_trace(wf,wf.mo_density_reconstruction.spin.packed);
}
void derive_molden_electron_counts(Wavefunction& wf) {
    if(wf.source!=WavefunctionSource::Molden)return;
    wf.alpha_electrons=wf.beta_electrons=0;wf.electron_counts_provenance=DataProvenance::Unavailable;
    const auto weights=occupation_weights(wf);
    wf.orbital_occupation_model=weights.total_info.occupation_model;
    wf.orbital_occupation_model_provenance=weights.total_info.model_provenance;
    if(weights.total_info.status!=NumericalStatus::Available || weights.spin_info.status!=NumericalStatus::Available)return;
    double alpha=0.0,beta=0.0;
    for(std::size_t i=0;i<weights.total.size();++i) {
        alpha+=(weights.total[i]+weights.spin[i])*0.5;beta+=(weights.total[i]-weights.spin[i])*0.5;
    }
    const auto integer=[](double value) {
        return std::isfinite(value) && value>=0.0 &&
            value<=static_cast<double>(std::numeric_limits<std::uint32_t>::max()) && std::floor(value)==value;
    };
    if(!integer(alpha) || !integer(beta))return;
    wf.alpha_electrons=static_cast<std::uint32_t>(alpha);wf.beta_electrons=static_cast<std::uint32_t>(beta);
    wf.electron_counts_provenance=DataProvenance::Derived;
}
std::vector<double> reconstruct_total_density_packed(const Wavefunction& wf) {
    return reconstruct_density(wf).total.packed;
}
std::vector<double> reconstruct_spin_density_packed(const Wavefunction& wf) {
    return reconstruct_density(wf).spin.packed;
}
} // namespace cov
