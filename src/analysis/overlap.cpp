#include "cov/overlap.hpp"
#include "cov/analysis_threads.hpp"
#include "cov/threading.hpp"

#include <algorithm>
#include <atomic>
#include <barrier>
#include <cmath>
#include <cstddef>
#include <latch>
#include <limits>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

#if defined(COV_ENABLE_CUDA)
#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#if defined(COV_CUDA_DYNAMIC_DENSE)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif
#endif

namespace cov {
namespace {

std::size_t dense_worker_count(const std::size_t n) noexcept {
    if (n<384u) return 1u;
    return analysis_thread_budget();
}

template<class Function>
void parallel_chunks(const std::size_t count,
                     const std::size_t workers,
                     Function&& function) {
    if (workers<=1u || count<workers) {
        function(0u,0u,count);
        return;
    }
    std::vector<cov::jthread> threads;
    threads.reserve(workers);
    for (std::size_t worker=0u;worker<workers;++worker) {
        threads.emplace_back([&,worker] {
            const std::size_t begin=count*worker/workers;
            const std::size_t end=count*(worker+1u)/workers;
            function(worker,begin,end);
        });
    }
}

double complete_metric_error_cpu(const std::vector<double>& c,
                                 const std::vector<double>& overlap,
                                 const std::size_t n,
                                 const std::size_t workers) {
    // Validate every entry of C^T*S*C, but keep only a small column block of
    // S*C live at once. The transposed coefficient block makes both dot
    // products contiguous and avoids the old large-matrix column sampling.
    std::vector<double> c_transposed(n*n,0.0);
    parallel_chunks(n,workers,[&](const std::size_t,
                                  const std::size_t begin,
                                  const std::size_t end) {
        for (std::size_t orbital=begin;orbital<end;++orbital) {
            double* output=c_transposed.data()+orbital*n;
            for (std::size_t basis=0u;basis<n;++basis) {
                output[basis]=c[basis*n+orbital];
            }
        }
    });

    constexpr std::size_t column_block=32u;
    double maximum_error=0.0;
    for (std::size_t column_begin=0u;column_begin<n;
         column_begin+=column_block) {
        const std::size_t block=std::min(column_block,n-column_begin);
        std::vector<double> sc_transposed(block*n,0.0);
        parallel_chunks(block*n,workers,[&](const std::size_t,
                                            const std::size_t begin,
                                            const std::size_t end) {
            for (std::size_t item=begin;item<end;++item) {
                const std::size_t local_column=item/n;
                const std::size_t row=item%n;
                const double* overlap_row=overlap.data()+row*n;
                const double* coefficient_column=
                    c_transposed.data()+(column_begin+local_column)*n;
                double value=0.0;
                for (std::size_t k=0u;k<n;++k) {
                    value+=overlap_row[k]*coefficient_column[k];
                }
                sc_transposed[local_column*n+row]=value;
            }
        });

        std::vector<double> worker_error(workers,0.0);
        parallel_chunks(n*block,workers,[&](const std::size_t worker,
                                            const std::size_t begin,
                                            const std::size_t end) {
            double local_error=0.0;
            for (std::size_t item=begin;item<end;++item) {
                const std::size_t orbital=item/block;
                const std::size_t local_column=item%block;
                const double* left=c_transposed.data()+orbital*n;
                const double* right=sc_transposed.data()+local_column*n;
                double value=0.0;
                for (std::size_t k=0u;k<n;++k) value+=left[k]*right[k];
                const std::size_t column=column_begin+local_column;
                const double target=orbital==column?1.0:0.0;
                local_error=std::max(
                    local_error,std::abs(value-target));
            }
            worker_error[worker]=local_error;
        });
        maximum_error=std::max(
            maximum_error,*std::max_element(
                worker_error.begin(),worker_error.end()));
    }
    return maximum_error;
}

#if defined(COV_ENABLE_CUDA)

class CudaDenseApi {
public:
    CudaDenseApi()=default;
    CudaDenseApi(const CudaDenseApi&)=delete;
    CudaDenseApi& operator=(const CudaDenseApi&)=delete;
    ~CudaDenseApi() {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        if (solver_module_) (void)FreeLibrary(solver_module_);
        if (blas_module_) (void)FreeLibrary(blas_module_);
#endif
    }

    bool load() {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        blas_module_=LoadLibraryW(L"cublas64_12.dll");
        solver_module_=LoadLibraryW(L"cusolver64_11.dll");
        if (!blas_module_ || !solver_module_) return false;
        cublas_create_=symbol<decltype(&cublasCreate_v2)>(
            blas_module_,"cublasCreate_v2");
        cublas_destroy_=symbol<decltype(&cublasDestroy_v2)>(
            blas_module_,"cublasDestroy_v2");
        cublas_dgemm_=symbol<decltype(&cublasDgemm_v2)>(
            blas_module_,"cublasDgemm_v2");
        solver_create_=symbol<decltype(&cusolverDnCreate)>(
            solver_module_,"cusolverDnCreate");
        solver_destroy_=symbol<decltype(&cusolverDnDestroy)>(
            solver_module_,"cusolverDnDestroy");
        solver_getrf_size_=symbol<decltype(&cusolverDnDgetrf_bufferSize)>(
            solver_module_,"cusolverDnDgetrf_bufferSize");
        solver_getrf_=symbol<decltype(&cusolverDnDgetrf)>(
            solver_module_,"cusolverDnDgetrf");
        solver_getrs_=symbol<decltype(&cusolverDnDgetrs)>(
            solver_module_,"cusolverDnDgetrs");
        return cublas_create_ && cublas_destroy_ && cublas_dgemm_ &&
               solver_create_ && solver_destroy_ && solver_getrf_size_ &&
               solver_getrf_ && solver_getrs_;
#else
        return true;
#endif
    }

    cublasStatus_t blas_create(cublasHandle_t* handle) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return cublas_create_(handle);
#else
        return cublasCreate(handle);
#endif
    }
    cublasStatus_t blas_destroy(cublasHandle_t handle) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return cublas_destroy_(handle);
#else
        return cublasDestroy(handle);
#endif
    }
    cublasStatus_t dgemm(cublasHandle_t handle,
                         cublasOperation_t transa,
                         cublasOperation_t transb,
                         int m,int n,int k,const double* alpha,
                         const double* a,int lda,const double* b,int ldb,
                         const double* beta,double* c,int ldc) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return cublas_dgemm_(handle,transa,transb,m,n,k,alpha,a,lda,b,ldb,
                             beta,c,ldc);
#else
        return cublasDgemm(handle,transa,transb,m,n,k,alpha,a,lda,b,ldb,
                           beta,c,ldc);
#endif
    }
    cusolverStatus_t solver_create(cusolverDnHandle_t* handle) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return solver_create_(handle);
#else
        return cusolverDnCreate(handle);
#endif
    }
    cusolverStatus_t solver_destroy(cusolverDnHandle_t handle) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return solver_destroy_(handle);
#else
        return cusolverDnDestroy(handle);
#endif
    }
    cusolverStatus_t getrf_buffer_size(cusolverDnHandle_t handle,int m,int n,
                                       double* a,int lda,int* size) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return solver_getrf_size_(handle,m,n,a,lda,size);
#else
        return cusolverDnDgetrf_bufferSize(handle,m,n,a,lda,size);
#endif
    }
    cusolverStatus_t getrf(cusolverDnHandle_t handle,int m,int n,double* a,
                           int lda,double* work,int* pivots,int* info) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return solver_getrf_(handle,m,n,a,lda,work,pivots,info);
#else
        return cusolverDnDgetrf(handle,m,n,a,lda,work,pivots,info);
#endif
    }
    cusolverStatus_t getrs(cusolverDnHandle_t handle,cublasOperation_t trans,
                           int n,int rhs,const double* a,int lda,
                           const int* pivots,double* b,int ldb,int* info) const {
#if defined(COV_CUDA_DYNAMIC_DENSE)
        return solver_getrs_(handle,trans,n,rhs,a,lda,pivots,b,ldb,info);
#else
        return cusolverDnDgetrs(
            handle,trans,n,rhs,a,lda,pivots,b,ldb,info);
#endif
    }

private:
#if defined(COV_CUDA_DYNAMIC_DENSE)
    template<class Function>
    static Function symbol(HMODULE module,const char* name) {
        return reinterpret_cast<Function>(GetProcAddress(module,name));
    }
    HMODULE blas_module_=nullptr;
    HMODULE solver_module_=nullptr;
    decltype(&cublasCreate_v2) cublas_create_=nullptr;
    decltype(&cublasDestroy_v2) cublas_destroy_=nullptr;
    decltype(&cublasDgemm_v2) cublas_dgemm_=nullptr;
    decltype(&cusolverDnCreate) solver_create_=nullptr;
    decltype(&cusolverDnDestroy) solver_destroy_=nullptr;
    decltype(&cusolverDnDgetrf_bufferSize) solver_getrf_size_=nullptr;
    decltype(&cusolverDnDgetrf) solver_getrf_=nullptr;
    decltype(&cusolverDnDgetrs) solver_getrs_=nullptr;
#endif
};

template<class T>
class DeviceBuffer {
public:
    DeviceBuffer()=default;
    DeviceBuffer(const DeviceBuffer&)=delete;
    DeviceBuffer& operator=(const DeviceBuffer&)=delete;
    ~DeviceBuffer() {
        if (data_) (void)cudaFree(data_);
    }
    bool allocate(const std::size_t count) {
        return cudaMalloc(reinterpret_cast<void**>(&data_),
                          count*sizeof(T))==cudaSuccess;
    }
    [[nodiscard]] T* get() const noexcept { return data_; }
private:
    T* data_=nullptr;
};

class CudaDenseHandles {
public:
    CudaDenseHandles()=default;
    CudaDenseHandles(const CudaDenseHandles&)=delete;
    CudaDenseHandles& operator=(const CudaDenseHandles&)=delete;
    ~CudaDenseHandles() {
        if (blas_) (void)api_.blas_destroy(blas_);
        if (solver_) (void)api_.solver_destroy(solver_);
    }
    bool create() {
        return api_.load() &&
               api_.solver_create(&solver_)==CUSOLVER_STATUS_SUCCESS &&
               api_.blas_create(&blas_)==CUBLAS_STATUS_SUCCESS;
    }
    [[nodiscard]] cusolverDnHandle_t solver() const noexcept {
        return solver_;
    }
    [[nodiscard]] cublasHandle_t blas() const noexcept { return blas_; }
    [[nodiscard]] const CudaDenseApi& api() const noexcept { return api_; }
private:
    CudaDenseApi api_;
    cusolverDnHandle_t solver_=nullptr;
    cublasHandle_t blas_=nullptr;
};

bool derive_overlap_cuda(const std::vector<double>& c,
                         const std::size_t n,
                         const double relative_tolerance,
                         std::vector<double>& overlap,
                         double& pivot_ratio,
                         double& maximum_error) {
    int device_count=0;
    if (cudaGetDeviceCount(&device_count)!=cudaSuccess || device_count<=0) {
        (void)cudaGetLastError();
        return false;
    }

    CudaDenseHandles handles;
    if (!handles.create()) return false;

    const std::size_t matrix_size=n*n;
    DeviceBuffer<double> matrix;
    DeviceBuffer<double> inverse_or_metric;
    DeviceBuffer<double> overlap_device;
    DeviceBuffer<double> temporary;
    DeviceBuffer<int> pivots;
    DeviceBuffer<int> info;
    if (!matrix.allocate(matrix_size) ||
        !inverse_or_metric.allocate(matrix_size) ||
        !overlap_device.allocate(matrix_size) ||
        !temporary.allocate(matrix_size) ||
        !pivots.allocate(n) || !info.allocate(1u)) {
        (void)cudaGetLastError();
        return false;
    }
    if (cudaMemcpy(matrix.get(),c.data(),matrix_size*sizeof(double),
                   cudaMemcpyHostToDevice)!=cudaSuccess) return false;

    int workspace_size=0;
    if (handles.api().getrf_buffer_size(
            handles.solver(),static_cast<int>(n),static_cast<int>(n),
            matrix.get(),static_cast<int>(n),&workspace_size)!=
        CUSOLVER_STATUS_SUCCESS || workspace_size<=0) return false;
    DeviceBuffer<double> workspace;
    if (!workspace.allocate(static_cast<std::size_t>(workspace_size))) {
        (void)cudaGetLastError();
        return false;
    }
    if (handles.api().getrf(
            handles.solver(),static_cast<int>(n),static_cast<int>(n),
            matrix.get(),static_cast<int>(n),workspace.get(),pivots.get(),
            info.get())!=CUSOLVER_STATUS_SUCCESS) return false;
    int host_info=0;
    if (cudaMemcpy(&host_info,info.get(),sizeof(int),
                   cudaMemcpyDeviceToHost)!=cudaSuccess || host_info!=0) {
        return false;
    }

    std::vector<double> lu(matrix_size,0.0);
    if (cudaMemcpy(lu.data(),matrix.get(),matrix_size*sizeof(double),
                   cudaMemcpyDeviceToHost)!=cudaSuccess) return false;
    double scale=0.0;
    for (const double value:c) scale=std::max(scale,std::abs(value));
    double minimum_pivot=std::numeric_limits<double>::infinity();
    double maximum_pivot=0.0;
    for (std::size_t i=0u;i<n;++i) {
        const double pivot=std::abs(lu[i*n+i]);
        if (!std::isfinite(pivot)) return false;
        minimum_pivot=std::min(minimum_pivot,pivot);
        maximum_pivot=std::max(maximum_pivot,pivot);
    }
    if (!(scale>0.0) || !std::isfinite(scale) ||
        minimum_pivot<=relative_tolerance*scale ||
        !(maximum_pivot>0.0)) return false;

    std::vector<double> identity(matrix_size,0.0);
    for (std::size_t i=0u;i<n;++i) identity[i*n+i]=1.0;
    if (cudaMemcpy(inverse_or_metric.get(),identity.data(),
                   matrix_size*sizeof(double),cudaMemcpyHostToDevice)!=
        cudaSuccess) return false;
    if (handles.api().getrs(
            handles.solver(),CUBLAS_OP_N,static_cast<int>(n),
            static_cast<int>(n),matrix.get(),static_cast<int>(n),pivots.get(),
            inverse_or_metric.get(),static_cast<int>(n),info.get())!=
        CUSOLVER_STATUS_SUCCESS) return false;
    if (cudaMemcpy(&host_info,info.get(),sizeof(int),
                   cudaMemcpyDeviceToHost)!=cudaSuccess || host_info!=0) {
        return false;
    }

    const double one=1.0;
    const double zero=0.0;
    // The uploaded row-major C is seen by cuBLAS as column-major C^T.
    // Solving (C^T)B=I gives B=C^-T, so B*B^T is C^-T*C^-1=S.
    if (handles.api().dgemm(
            handles.blas(),CUBLAS_OP_N,CUBLAS_OP_T,
            static_cast<int>(n),static_cast<int>(n),static_cast<int>(n),
            &one,inverse_or_metric.get(),static_cast<int>(n),
            inverse_or_metric.get(),static_cast<int>(n),&zero,
            overlap_device.get(),static_cast<int>(n))!=
        CUBLAS_STATUS_SUCCESS) return false;

    // Restore C^T and form the complete C^T*S*C matrix with two GEMMs. No
    // orbital-column sampling is allowed: the reported error is the true
    // maximum over all n^2 entries.
    if (cudaMemcpy(matrix.get(),c.data(),matrix_size*sizeof(double),
                   cudaMemcpyHostToDevice)!=cudaSuccess) return false;
    if (handles.api().dgemm(
            handles.blas(),CUBLAS_OP_N,CUBLAS_OP_T,
            static_cast<int>(n),static_cast<int>(n),static_cast<int>(n),
            &one,overlap_device.get(),static_cast<int>(n),
            matrix.get(),static_cast<int>(n),&zero,temporary.get(),
            static_cast<int>(n))!=CUBLAS_STATUS_SUCCESS) return false;
    if (handles.api().dgemm(
            handles.blas(),CUBLAS_OP_N,CUBLAS_OP_N,
            static_cast<int>(n),static_cast<int>(n),static_cast<int>(n),
            &one,matrix.get(),static_cast<int>(n),temporary.get(),
            static_cast<int>(n),&zero,inverse_or_metric.get(),
            static_cast<int>(n))!=CUBLAS_STATUS_SUCCESS) return false;

    std::vector<double> metric(matrix_size,0.0);
    overlap.assign(matrix_size,0.0);
    if (cudaMemcpy(overlap.data(),overlap_device.get(),
                   matrix_size*sizeof(double),cudaMemcpyDeviceToHost)!=
            cudaSuccess ||
        cudaMemcpy(metric.data(),inverse_or_metric.get(),
                   matrix_size*sizeof(double),cudaMemcpyDeviceToHost)!=
            cudaSuccess) {
        overlap.clear();
        return false;
    }
    maximum_error=0.0;
    for (std::size_t row=0u;row<n;++row) {
        for (std::size_t column=0u;column<n;++column) {
            const double target=row==column?1.0:0.0;
            maximum_error=std::max(
                maximum_error,
                std::abs(metric[column*n+row]-target));
        }
    }
    pivot_ratio=minimum_pivot/maximum_pivot;
    return std::isfinite(maximum_error);
}

#endif

std::vector<const MolecularOrbital*> complete_spin_block(const Wavefunction& wf,
                                                         const Spin spin) {
    std::vector<const MolecularOrbital*> block;
    block.reserve(wf.basis_count);
    for (const auto& mo : wf.orbitals) {
        if (mo.spin == spin) block.push_back(&mo);
    }
    if (block.size() != wf.basis_count) return {};
    for (const auto* mo : block) {
        if (!mo || mo->coefficients.size() != wf.basis_count) return {};
    }
    return block;
}

bool invert_square(std::vector<double>& a,
                   std::vector<double>& inverse,
                   const std::size_t n,
                   const double relative_tolerance,
                   double& pivot_ratio) {
    inverse.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) inverse[i * n + i] = 1.0;

    double scale = 0.0;
    for (const double value : a) scale = std::max(scale, std::abs(value));
    if (!(scale > 0.0) || !std::isfinite(scale)) return false;

    double min_pivot = std::numeric_limits<double>::infinity();
    double max_pivot = 0.0;

    const std::size_t workers=dense_worker_count(n);
    std::atomic<bool> stop{false};
    std::size_t active_column=0u;
    std::vector<cov::jthread> elimination_threads;
    std::unique_ptr<std::barrier<>> phase;
    std::latch launch_gate(1u);
    if (workers>1u) {
        try {
            phase=std::make_unique<std::barrier<>>(
                static_cast<std::ptrdiff_t>(workers+1u));
            elimination_threads.reserve(workers);
            for (std::size_t worker=0u;worker<workers;++worker) {
                elimination_threads.emplace_back([&,worker] {
                    // Do not enter a fixed-participant barrier until every
                    // worker has been constructed. If construction throws,
                    // the catch path opens this gate and all partial workers
                    // leave without waiting for participants that never exist.
                    launch_gate.wait();
                    if (stop.load(std::memory_order_acquire)) return;
                    for (;;) {
                        phase->arrive_and_wait();
                        if (stop.load(std::memory_order_acquire)) break;
                        const std::size_t col=active_column;
                        const std::size_t begin=n*worker/workers;
                        const std::size_t end=n*(worker+1u)/workers;
                        for (std::size_t row=begin;row<end;++row) {
                            if (row==col) continue;
                            const double factor=a[row*n+col];
                            if (factor==0.0) continue;
                            a[row*n+col]=0.0;
                            for (std::size_t j=col+1u;j<n;++j) {
                                a[row*n+j]-=factor*a[col*n+j];
                            }
                            for (std::size_t j=0u;j<n;++j) {
                                inverse[row*n+j]-=
                                    factor*inverse[col*n+j];
                            }
                        }
                        phase->arrive_and_wait();
                    }
                });
            }
        } catch (...) {
            stop.store(true,std::memory_order_release);
            launch_gate.count_down();
            elimination_threads.clear();
            return false;
        }
        launch_gate.count_down();
    }

    bool valid=true;
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        double pivot_abs = std::abs(a[col * n + col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const double candidate = std::abs(a[row * n + col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if (!std::isfinite(pivot_abs) ||
            pivot_abs <= relative_tolerance * scale) {
            valid=false;
            break;
        }

        min_pivot = std::min(min_pivot, pivot_abs);
        max_pivot = std::max(max_pivot, pivot_abs);
        if (pivot_row != col) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(a[col * n + j], a[pivot_row * n + j]);
                std::swap(inverse[col * n + j], inverse[pivot_row * n + j]);
            }
        }

        const double pivot = a[col * n + col];
        for (std::size_t j = 0; j < n; ++j) {
            a[col * n + j] /= pivot;
            inverse[col * n + j] /= pivot;
        }

        if (workers>1u) {
            active_column=col;
            phase->arrive_and_wait();
            phase->arrive_and_wait();
        } else {
            for (std::size_t row = 0; row < n; ++row) {
                if (row == col) continue;
                const double factor = a[row * n + col];
                if (factor == 0.0) continue;
                a[row * n + col] = 0.0;
                for (std::size_t j = col + 1; j < n; ++j) {
                    a[row * n + j] -= factor * a[col * n + j];
                }
                for (std::size_t j = 0; j < n; ++j) {
                    inverse[row * n + j] -= factor * inverse[col * n + j];
                }
            }
        }
    }

    if (workers>1u) {
        stop.store(true,std::memory_order_release);
        phase->arrive_and_wait();
        elimination_threads.clear();
    }
    if (!valid) return false;

    pivot_ratio = max_pivot > 0.0 ? min_pivot / max_pivot : 0.0;
    return true;
}

} // namespace

OverlapDerivationResult derive_ao_overlap_from_mos(const Wavefunction& wavefunction,
                                                   const double relative_pivot_tolerance,
                                                   const std::size_t maximum_basis) {
    OverlapDerivationResult result;
    const std::size_t n = wavefunction.basis_count;
    if (n == 0 || n > maximum_basis) return result;

    auto block = complete_spin_block(wavefunction, Spin::Alpha);
    if (block.empty()) block = complete_spin_block(wavefunction, Spin::Beta);
    if (block.empty()) return result;

    // C is stored AO-major here: C(mu,i).
    std::vector<double> c(n * n, 0.0);
    for (std::size_t orbital = 0; orbital < n; ++orbital) {
        for (std::size_t basis = 0; basis < n; ++basis) {
            c[basis * n + orbital] =
                static_cast<double>(block[orbital]->coefficients[basis]);
        }
    }

#if defined(COV_ENABLE_CUDA)
    // Dense recovery becomes GPU-efficient at this scale. Every CUDA API
    // failure, a missing driver/device, allocation pressure, or a numerically
    // rejected LU factorization falls through to the portable CPU path.
    if (n>=384u) {
        std::vector<double> cuda_overlap;
        double cuda_pivot_ratio=0.0;
        double cuda_error=0.0;
        if (derive_overlap_cuda(c,n,relative_pivot_tolerance,cuda_overlap,
                                cuda_pivot_ratio,cuda_error) &&
            cuda_error<=5.0e-6) {
            result.matrix=std::move(cuda_overlap);
            result.pivot_ratio=cuda_pivot_ratio;
            result.max_orthonormality_error=cuda_error;
            return result;
        }
    }
#endif

    std::vector<double> work = c;
    std::vector<double> c_inverse;
    if (!invert_square(work, c_inverse, n, relative_pivot_tolerance,
                       result.pivot_ratio)) {
        result.pivot_ratio = 0.0;
        return result;
    }

    const std::size_t workers=dense_worker_count(n);

    // If X=C^{-1}, then S=X^T X. Accumulate one output row per worker with
    // the inner loop contiguous in both X and S; the previous column-dot
    // ordering was strided and made 900+ basis recovery look like a hang.
    result.matrix.assign(n * n, 0.0);
    parallel_chunks(n,workers,[&](const std::size_t,
                                  const std::size_t begin,
                                  const std::size_t end) {
        for (std::size_t mu=begin;mu<end;++mu) {
            double* output=result.matrix.data()+mu*n;
            for (std::size_t i=0u;i<n;++i) {
                const double scale=c_inverse[i*n+mu];
                const double* inverse_row=c_inverse.data()+i*n;
                for (std::size_t nu=0u;nu<=mu;++nu) {
                    output[nu]+=scale*inverse_row[nu];
                }
            }
        }
    });
    parallel_chunks(n,workers,[&](const std::size_t,
                                  const std::size_t begin,
                                  const std::size_t end) {
        for (std::size_t mu=begin;mu<end;++mu) {
            for (std::size_t nu=0u;nu<mu;++nu) {
                result.matrix[nu*n+mu]=result.matrix[mu*n+nu];
            }
        }
    });

    const double max_error=complete_metric_error_cpu(
        c,result.matrix,n,workers);
    result.max_orthonormality_error = max_error;

    if (!std::isfinite(max_error) || max_error > 5.0e-6) {
        result.matrix.clear();
    }
    return result;
}

} // namespace cov
