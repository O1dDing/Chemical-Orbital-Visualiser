#include "cov/molecular_point_group_frame.hpp"
#include "cov/nbo_salc.hpp"
#include <Eigen/Core>
#include "cov/orbital_symmetry.hpp"

// Keep the already-tested finite/linear machinery in one translation unit so
// the final dispatcher can reuse its AO-operation helpers without exposing them
// as public API. The included implementation is renamed here; CMake compiles
// this file instead of compiling orbital_symmetry_dispatch.cpp directly.
#define derive_orbital_symmetry derive_orbital_symmetry_dispatch_v1
#include "orbital_symmetry_dispatch.cpp"
#undef derive_orbital_symmetry

namespace cov {
OrbitalSymmetryResult derive_orbital_symmetry(
    Wavefunction& wavefunction,
    const OrbitalSymmetryOptions& options) {
    if (wavefunction.orbitals.empty() || wavefunction.basis_count == 0 ||
        wavefunction.ao_overlap.size() !=
            static_cast<std::size_t>(wavefunction.basis_count) *
                wavefunction.basis_count) {
        return {};
    }

    const MolecularSymmetry symmetry = analyse_molecular_symmetry(wavefunction);

    if(!symmetry.linear){
        OrbitalSymmetryResult result;result.point_group=symmetry.point_group;
        const auto complete=complete_molecular_point_group(wavefunction,symmetry);
        const auto table=molecular_point_group_irreps(wavefunction,complete);
        if(!table.valid)return result;
        using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
        const auto n=std::size_t(wavefunction.basis_count),m=wavefunction.orbitals.size();Matrix q=Matrix::Zero(n,m);
        for(std::size_t c=0;c<m;++c){if(wavefunction.orbitals[c].coefficients.size()!=n)return result;
            for(std::size_t r=0;r<n;++r)q(r,c)=wavefunction.orbitals[c].coefficients[r];}
        const Matrix sq=Eigen::Map<const Matrix>(wavefunction.ao_overlap.data(),n,n)*q;
        const auto groups=energy_groups(wavefunction,options.degeneracy_tolerance_hartree);
        std::vector<std::vector<double>> measured(groups.size());std::vector<double> retention(groups.size(),1);
        std::vector<bool> eligible(groups.size());for(std::size_t b=0;b<groups.size();++b){eligible[b]=group_unlabelled(wavefunction,groups[b]);if(eligible[b])++result.groups_examined;}
        const std::vector<double> packed(q.data(),q.data()+q.size());
        for(const auto& operation:complete.operations){const auto transformed=apply_orbital_symmetry_operation(wavefunction,operation,packed,m);
            if(transformed.size()!=n*m)return result;const Eigen::Map<const Matrix> tq(transformed.data(),n,m);
            const Matrix d=sq.transpose()*tq;const Matrix stq=Eigen::Map<const Matrix>(wavefunction.ao_overlap.data(),n,n)*tq;
            for(std::size_t b=0;b<groups.size();++b){if(!eligible[b])continue;double chi=0;
                for(auto i:groups[b]){chi+=d(i,i);double projected=0;for(auto j:groups[b])projected+=d(j,i)*d(j,i);
                    const double norm=tq.col(i).dot(stq.col(i));if(!std::isfinite(norm)||norm<=1e-12){eligible[b]=false;break;}
                    retention[b]=std::min(retention[b],projected/norm);}
                measured[b].push_back(chi);
            }
        }
        for(std::size_t b=0;b<groups.size();++b){if(!eligible[b])continue;result.worst_subspace_retention=std::min(result.worst_subspace_retention,retention[b]);
            if(retention[b]<options.minimum_subspace_retention)continue;
            std::string label;double maximum_error=0;for(const auto& row:table.rows){if(row.dimension!=groups[b].size()||row.characters.size()!=measured[b].size())continue;
                bool same=true;for(std::size_t g=0;g<row.characters.size();++g)if(!std::isfinite(measured[b][g])||std::abs(measured[b][g]-row.characters[g])>options.character_tolerance){same=false;break;}
                if(same){if(!label.empty()){label.clear();break;}label=row.label;for(std::size_t g=0;g<row.characters.size();++g)maximum_error=std::max(maximum_error,std::abs(measured[b][g]-row.characters[g]));}}
            if(label.empty())continue;
            DerivedOrbitalSymmetryAssignment evidence;evidence.point_group=symmetry.point_group;evidence.label=label;evidence.orbital_indices=groups[b];evidence.subspace_retention=retention[b];evidence.centre_bohr=symmetry.centre_bohr;
            evidence.maximum_character_error=maximum_error;evidence.axes_available=true;evidence.principal_axis=table.principal_axis;evidence.secondary_axis=table.reference_axis;evidence.axis_convention=table.axis_detail;
            wavefunction.derived_orbital_symmetry_assignments.push_back(std::move(evidence));++result.groups_labelled;
            for(auto i:groups[b]){wavefunction.orbitals[i].symmetry=label;wavefunction.orbitals[i].symmetry_provenance=DataProvenance::Derived;++result.orbitals_labelled;}
        }
        return result;
    }

    return derive_orbital_symmetry_dispatch_v1(wavefunction, options);
}

} // namespace cov
