#include "cov/overlap.hpp"
#include "cov/ao_angular_basis.hpp"
#include "cov/analysis_threads.hpp"
#include "cov/threading.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <exception>
#include <map>
#include <mutex>
#include <numbers>
#include <stdexcept>
#include <thread>
#include <utility>

namespace cov {
namespace {
using ao_angular::Powers;
using ao_angular::Term;
using ao_angular::Polynomial;
using ao_angular::cartesian_components;
using ao_angular::solid_harmonic;
using ao_angular::odd_factorial;
struct PreparedPrimitive {double exponent=0.0;double weight=0.0;};
struct PreparedShell {
    std::array<double,3> centre{};
    std::size_t offset=0;
    int l=0;
    std::vector<Polynomial> components;
    std::vector<PreparedPrimitive> primitives;
};

std::vector<PreparedShell> prepare(const Wavefunction& wf) {
    std::vector<PreparedShell> result;
    std::vector<bool> covered(wf.basis_count,false);
    for (const auto& shell:wf.shells) {
        if (shell.atom_index>=wf.atoms.size() || shell.angular_momentum>4 ||
            shell.primitive_count==0 || shell.primitive_offset>wf.primitives.size() ||
            shell.primitive_count>wf.primitives.size()-shell.primitive_offset) {
            throw std::invalid_argument("Invalid AO shell/primitive/atom range");
        }
        const auto count=static_cast<std::size_t>(shell_basis_count(shell));
        if (shell.basis_offset>wf.basis_count || count>wf.basis_count-shell.basis_offset)
            throw std::invalid_argument("AO shell range exceeds the declared basis");
        for (std::size_t i=0;i<count;++i) {
            const auto index=shell.basis_offset+i;
            if (covered[index]) throw std::invalid_argument("AO shell ranges overlap");
            covered[index]=true;
        }
        const auto& atom=wf.atoms[shell.atom_index];
        PreparedShell item;
        item.centre={atom.x,atom.y,atom.z};item.offset=shell.basis_offset;item.l=shell.angular_momentum;
        for (const auto x:item.centre) if (!std::isfinite(x)) throw std::invalid_argument("Nonfinite AO centre");
        if (shell.pure) {
            for (std::size_t i=0;i<count;++i) item.components.push_back(solid_harmonic(item.l,static_cast<int>(i)));
        } else {
            for (const auto& powers:cartesian_components(item.l)) {
                const double factor=1.0/std::sqrt(odd_factorial(2*powers[0]-1)*
                    odd_factorial(2*powers[1]-1)*odd_factorial(2*powers[2]-1));
                item.components.push_back({{powers,factor}});
            }
        }
        for (std::size_t i=0;i<shell.primitive_count;++i) {
            const auto& primitive=wf.primitives[shell.primitive_offset+i];
            const double alpha=primitive.exponent;
            if (!std::isfinite(alpha) || alpha<=0 || !std::isfinite(primitive.coefficient))
                throw std::invalid_argument("Invalid primitive exponent or contraction coefficient");
            const double normalization=shell.pure
                ? std::sqrt(2.0*std::pow(2.0*alpha,item.l+1.5)/std::tgamma(item.l+1.5))
                : std::pow(2.0*alpha/std::numbers::pi,0.75)*std::pow(4.0*alpha,item.l*0.5);
            const double weight=primitive.coefficient*normalization;
            if (!std::isfinite(weight)) throw std::overflow_error("Primitive normalization exceeds numeric range");
            item.primitives.push_back({alpha,weight});
        }
        result.push_back(std::move(item));
    }
    if (std::find(covered.begin(),covered.end(),false)!=covered.end())
        throw std::invalid_argument("Declared AO basis contains an unrepresented direction");
    return result;
}

using Moments=std::array<std::array<double,5>,5>;
Moments moments(const int left,const int right,const double pa,const double pb,const double p) {
    Moments result{};
    result[0][0]=std::sqrt(std::numbers::pi/p);
    const double reciprocal=0.5/p;
    for (int i=1;i<=left;++i)
        result[i][0]=pa*result[i-1][0]+(i>1?(i-1)*reciprocal*result[i-2][0]:0.0);
    for (int j=1;j<=right;++j) {
        for (int i=0;i<=left;++i) {
            result[i][j]=pb*result[i][j-1]+(j>1?(j-1)*reciprocal*result[i][j-2]:0.0)+
                (i>0?i*reciprocal*result[i-1][j-1]:0.0);
        }
    }
    return result;
}

void integrate_pair(const PreparedShell& left,const PreparedShell& right,
                    std::vector<double>& matrix,const std::size_t dimension,const bool same_shell) {
    const auto nl=left.components.size(),nr=right.components.size();
    std::vector<double> block(nl*nr,0.0);
    double distance2=0.0;
    for (int axis=0;axis<3;++axis) {
        const double difference=left.centre[axis]-right.centre[axis];distance2+=difference*difference;
    }
    for (const auto& a:left.primitives) {
        for (const auto& b:right.primitives) {
            const double p=a.exponent+b.exponent;
            const double factor=a.weight*b.weight*std::exp(-a.exponent*(b.exponent/p)*distance2);
            std::array<Moments,3> integral;
            for (int axis=0;axis<3;++axis) {
                // Differences avoid subtracting large translated centres.
                const double displacement=right.centre[axis]-left.centre[axis];
                integral[axis]=moments(left.l,right.l,(b.exponent/p)*displacement,
                                      -(a.exponent/p)*displacement,p);
            }
            for (std::size_t i=0;i<nl;++i) {
                const auto end=same_shell?i+1:nr;
                for (std::size_t j=0;j<end;++j) {
                    double angular=0.0;
                    for (const auto& x:left.components[i]) {
                        for (const auto& y:right.components[j]) {
                            double value=x.coefficient*y.coefficient;
                            for (int axis=0;axis<3;++axis) value*=integral[axis][x.powers[axis]][y.powers[axis]];
                            angular+=value;
                        }
                    }
                    block[i*nr+j]+=factor*angular;
                }
            }
        }
    }
    for (std::size_t i=0;i<nl;++i) {
        const auto end=same_shell?i+1:nr;
        for (std::size_t j=0;j<end;++j) {
            matrix[(left.offset+i)*dimension+right.offset+j]=block[i*nr+j];
            matrix[(right.offset+j)*dimension+left.offset+i]=block[i*nr+j];
        }
    }
}
}

BasisOverlapResult derive_ao_overlap_from_basis(const Wavefunction& wf) {
    BasisOverlapResult result;
    if (wf.basis_count==0 || wf.shells.empty() || wf.atoms.empty() || wf.primitives.empty()) {
        result.status=NumericalStatus::MissingInput;result.detail="AO basis data are absent";return result;
    }
    try {
        const auto shells=prepare(wf);
        const auto n=static_cast<std::size_t>(wf.basis_count);
        if (n>result.matrix.max_size()/n) throw std::length_error("AO matrix size exceeds addressable storage");
        result.matrix.assign(n*n,0.0);
        std::vector<std::pair<std::size_t,std::size_t>> pairs;
        for (std::size_t i=0;i<shells.size();++i) for (std::size_t j=0;j<=i;++j) pairs.emplace_back(i,j);
        std::atomic<std::size_t> next{0};
        std::atomic<bool> failed{false};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        const auto work=[&] {
            try {
                while (!failed.load()) {
                    const auto index=next.fetch_add(1);
                    if (index>=pairs.size()) return;
                    const auto [i,j]=pairs[index];integrate_pair(shells[i],shells[j],result.matrix,n,i==j);
                }
            } catch (...) {
                failed.store(true);
                const std::lock_guard lock(error_mutex);
                if (!worker_error) worker_error=std::current_exception();
            }
        };
        const auto workers=std::min(analysis_thread_budget(),pairs.size());
        {
            std::vector<cov::jthread> threads;
            for (std::size_t i=1;i<workers;++i) threads.emplace_back(work);
            work();
        }
        if (worker_error) std::rethrow_exception(worker_error);
        if (!std::all_of(result.matrix.begin(),result.matrix.end(),[](double x){return std::isfinite(x);}))
            throw std::overflow_error("Nonfinite analytic AO integral");
        result.status=NumericalStatus::Available;
        result.detail="Contracted normalized Gaussian basis integrals; COV Cartesian/Condon-Shortley real spherical convention; no MO fit";
    } catch (const std::invalid_argument& error) {
        result.matrix.clear();result.status=NumericalStatus::InvalidInput;result.detail=error.what();
    } catch (const std::exception& error) {
        result.matrix.clear();result.status=NumericalStatus::Failed;result.detail=error.what();
    }
    return result;
}

AoMetricDiagnostics inspect_ao_metric(const Wavefunction& wf,const std::vector<double>& matrix) {
    AoMetricDiagnostics result;
    const auto n=static_cast<std::size_t>(wf.basis_count);result.dimension=n;
    if (n==0 || matrix.empty()) {result.status=NumericalStatus::MissingInput;result.detail="AO matrix absent";return result;}
    if (matrix.size()/n!=n || matrix.size()%n!=0) {
        result.status=NumericalStatus::InvalidInput;result.detail="AO matrix dimension mismatch";return result;
    }
    result.finite=std::all_of(matrix.begin(),matrix.end(),[](double x){return std::isfinite(x);});
    if (!result.finite) {result.status=NumericalStatus::InvalidInput;result.detail="Nonfinite AO matrix";return result;}
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    const Eigen::Map<const Matrix> metric(matrix.data(),static_cast<Eigen::Index>(n),static_cast<Eigen::Index>(n));
    const double scale=std::max(1.0,metric.cwiseAbs().maxCoeff());
    result.symmetry_error=(metric-metric.transpose()).cwiseAbs().maxCoeff();
    const double roundoff=64.0*std::numeric_limits<double>::epsilon()*static_cast<double>(n)*scale;
    if (result.symmetry_error>roundoff) {result.status=NumericalStatus::InvalidInput;result.detail="Nonsymmetric AO matrix";return result;}
    if ((metric.diagonal().array()<=0.0).any()) {
        result.status=NumericalStatus::InvalidInput;result.detail="AO metric contains a zero or negative basis norm";return result;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(metric);
    if (solver.info()!=Eigen::Success) {result.status=NumericalStatus::Failed;result.detail="AO eigensolver did not converge";return result;}
    const auto& values=solver.eigenvalues();
    result.minimum_eigenvalue=values[0];result.maximum_eigenvalue=values[values.size()-1];
    result.numerical_rank_tolerance=64.0*std::numeric_limits<double>::epsilon()*static_cast<double>(n)*
        std::max(1.0,std::abs(result.maximum_eigenvalue));
    result.numerical_rank=static_cast<std::size_t>((values.array()>result.numerical_rank_tolerance).count());
    result.positive_semidefinite=result.minimum_eigenvalue>=-result.numerical_rank_tolerance;
    if (result.minimum_eigenvalue>0) result.condition_number=result.maximum_eigenvalue/result.minimum_eigenvalue;
    result.eigen_residual=(metric*solver.eigenvectors()-solver.eigenvectors()*values.asDiagonal()).cwiseAbs().maxCoeff()/scale;
    if (!std::isfinite(result.eigen_residual) || result.eigen_residual>64.0*std::numeric_limits<double>::epsilon()*static_cast<double>(n)) {
        result.status=NumericalStatus::Failed;result.detail="AO spectral decomposition residual failed";return result;
    }
    result.status=result.positive_semidefinite?NumericalStatus::Available:NumericalStatus::InvalidInput;
    result.detail=result.positive_semidefinite?"Full AO spectrum; numerical rank reported without truncation":"Indefinite AO metric";
    for (const auto spin:{Spin::Alpha,Spin::Beta}) {
        std::vector<const MolecularOrbital*> orbitals;
        for (const auto& mo:wf.orbitals) if (mo.spin==spin) orbitals.push_back(&mo);
        OrbitalMetricDiagnostics block;block.spin=spin;block.orbital_count=orbitals.size();
        if (orbitals.empty()) {
            block.status=NumericalStatus::MissingInput;block.detail="No explicit MO block for this spin";
            result.orbital_blocks.push_back(block);continue;
        }
        const auto m=static_cast<Eigen::Index>(orbitals.size());
        Eigen::MatrixXd coefficients(static_cast<Eigen::Index>(n),m);
        bool valid=true;
        for (Eigen::Index j=0;j<m;++j) {
            const auto& column=orbitals[static_cast<std::size_t>(j)]->coefficients;
            if (column.size()!=n || !std::all_of(column.begin(),column.end(),[](double x){return std::isfinite(x);})) {valid=false;break;}
            for (std::size_t i=0;i<n;++i) coefficients(static_cast<Eigen::Index>(i),j)=column[i];
        }
        if (!valid) {
            block.status=NumericalStatus::InvalidInput;block.detail="MO coefficient dimension mismatch or nonfinite coefficient";
            result.orbital_blocks.push_back(block);continue;
        }
        const Eigen::MatrixXd gram=coefficients.transpose()*metric*coefficients;
        if (!gram.allFinite()) {
            block.status=NumericalStatus::Failed;block.detail="Nonfinite MO metric product";
            result.orbital_blocks.push_back(block);continue;
        }
        block.maximum_diagonal_error=(gram.diagonal().array()-1.0).abs().maxCoeff();
        if (m>1) {
            block.maximum_off_diagonal_error=0.0;
            for (Eigen::Index i=0;i<m;++i) for (Eigen::Index j=0;j<i;++j)
                block.maximum_off_diagonal_error=std::max(block.maximum_off_diagonal_error,std::abs(gram(i,j)));
        }
        const bool orthonormal=block.maximum_diagonal_error<=block.orthonormality_tolerance &&
            (m==1 || block.maximum_off_diagonal_error<=block.orthonormality_tolerance);
        block.status=orthonormal?NumericalStatus::Available:NumericalStatus::InvalidInput;
        block.detail=orthonormal?"All entries of the spin MO Gram matrix satisfy the recorded tolerance":
            "Input MO block is not orthonormal in this AO metric; coefficients retained unchanged";
        result.orbital_blocks.push_back(block);
    }
    return result;
}

void establish_ao_metric(Wavefunction& wf) {
    // Preserve the producer matrix even when it is invalid, but never use it
    // as authority over integrals of the actual supplied Gaussian functions.
    if (wf.producer_ao_overlap.empty() && wf.ao_overlap_provenance==DataProvenance::Producer)
        wf.producer_ao_overlap=wf.ao_overlap;
    wf.producer_ao_metric_diagnostics=inspect_ao_metric(wf,wf.producer_ao_overlap);
    wf.ao_overlap.clear();
    wf.ao_overlap_provenance=DataProvenance::Unavailable;
    wf.ao_overlap_orthonormality_error=std::numeric_limits<double>::quiet_NaN();
    wf.producer_ao_overlap_basis_error=std::numeric_limits<double>::quiet_NaN();
    wf.producer_ao_overlap_basis_status=wf.producer_ao_overlap.empty()
        ?NumericalStatus::MissingInput:NumericalStatus::NotComputed;
    auto basis=derive_ao_overlap_from_basis(wf);
    if (!basis.available()) {
        wf.ao_metric_diagnostics={};
        wf.ao_metric_diagnostics.status=basis.status;
        wf.ao_metric_diagnostics.dimension=wf.basis_count;
        wf.ao_metric_diagnostics.detail=basis.detail;
        return;
    }
    wf.ao_overlap=std::move(basis.matrix);
    wf.ao_overlap_provenance=DataProvenance::Derived;
    wf.ao_metric_diagnostics=inspect_ao_metric(wf,wf.ao_overlap);
    wf.ao_metric_diagnostics.detail=basis.detail+"; "+wf.ao_metric_diagnostics.detail;
    bool measured=false;
    double error=0.0;
    for (const auto& block:wf.ao_metric_diagnostics.orbital_blocks) {
        if (std::isfinite(block.maximum_diagonal_error)) {
            measured=true;error=std::max(error,block.maximum_diagonal_error);
        }
        if (std::isfinite(block.maximum_off_diagonal_error))
            error=std::max(error,block.maximum_off_diagonal_error);
    }
    if (measured) wf.ao_overlap_orthonormality_error=error;
    if (!wf.producer_ao_overlap.empty()) {
        if (wf.producer_ao_overlap.size()!=wf.ao_overlap.size() || !wf.producer_ao_metric_diagnostics.finite) {
            wf.producer_ao_overlap_basis_status=NumericalStatus::InvalidInput;
        } else {
            double difference=0.0;
            for (std::size_t i=0;i<wf.ao_overlap.size();++i)
                difference=std::max(difference,std::abs(wf.ao_overlap[i]-wf.producer_ao_overlap[i]));
            wf.producer_ao_overlap_basis_error=difference;
            wf.producer_ao_overlap_basis_status=difference<=wf.producer_ao_overlap_basis_tolerance
                ?NumericalStatus::Available:NumericalStatus::InvalidInput;
        }
    }
}

bool orbital_metric_usable(const Wavefunction& wf) noexcept {
    if (wf.ao_metric_diagnostics.status!=NumericalStatus::Available || wf.orbitals.empty()) return false;
    bool present=false;
    for (const auto& block:wf.ao_metric_diagnostics.orbital_blocks) {
        if (block.orbital_count==0) continue;
        present=true;
        if (block.status!=NumericalStatus::Available) return false;
    }
    return present;
}

} // namespace cov
