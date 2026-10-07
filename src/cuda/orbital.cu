#include "cov/cuda_orbital.hpp"
#include "cov/gl_api.hpp"

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif
#include <GL/gl.h>

#include <cuda_gl_interop.h>
#include <cuda_runtime.h>

#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef GL_TEXTURE_3D
#define GL_TEXTURE_3D 0x806F
#endif
#ifndef GL_TEXTURE_BINDING_3D
#define GL_TEXTURE_BINDING_3D 0x806A
#endif
#ifndef GL_RED
#define GL_RED 0x1903
#endif

namespace cov {
namespace {

constexpr float kPi = 3.14159265358979323846f;

struct GpuPrimitive {
    float exponent;
    float coefficient;
};

struct GpuShell {
    float cx;
    float cy;
    float cz;
    std::uint32_t primitive_offset;
    std::uint32_t primitive_count;
    std::uint32_t basis_offset;
    std::uint8_t l;
    std::uint8_t pure;
    std::uint16_t pad;
};

float pack_gpu_scalar(const double value, const char* field) {
    const float packed=static_cast<float>(value);
    if (!std::isfinite(value) || !std::isfinite(packed)) {
        throw std::runtime_error(std::string("Nonfinite or overflowing GPU input: ")+field);
    }
    return packed;
}

void cuda_check(const cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(status));
    }
}

__device__ __forceinline__ float powi(float x, int n) {
    float out = 1.0f;
    for (int i = 0; i < n; ++i) out *= x;
    return out;
}

__device__ __forceinline__ int double_factorial_odd(int n) {
    if (n <= 0) return 1;
    int v = 1;
    for (int k = n; k > 1; k -= 2) v *= k;
    return v;
}

__device__ __forceinline__ float cartesian_primitive_norm(
    const float alpha, const int ax, const int ay, const int az) {
    const int l = ax + ay + az;
    const float base = powf(2.0f * alpha / kPi, 0.75f);
    const float angular = powf(4.0f * alpha, 0.5f * static_cast<float>(l));
    const float denom = sqrtf(static_cast<float>(
        double_factorial_odd(2 * ax - 1) *
        double_factorial_odd(2 * ay - 1) *
        double_factorial_odd(2 * az - 1)));
    return base * angular / denom;
}

__device__ __forceinline__ float gamma_l_plus_three_halves(const int l) {
    const float sqrt_pi = sqrtf(kPi);
    switch (l) {
        case 0: return 0.5f * sqrt_pi;
        case 1: return 0.75f * sqrt_pi;
        case 2: return 1.875f * sqrt_pi;
        case 3: return 6.5625f * sqrt_pi;
        default: return 29.53125f * sqrt_pi;
    }
}

__device__ __forceinline__ float spherical_primitive_norm(const float alpha, const int l) {
    const float power = powf(2.0f * alpha, static_cast<float>(l) + 1.5f);
    return sqrtf(2.0f * power / gamma_l_plus_three_halves(l));
}

__device__ __forceinline__ float real_solid_harmonic(
    const int l, const int index, const float x, const float y, const float z) {
    const float x2 = x * x;
    const float y2 = y * y;
    const float z2 = z * z;
    const float r2 = x * x + y * y + z * z;
    if (l == 0) return 0.28209479177387814f;

    // Normalised real solid harmonics in Molden order:
    // m = 0, +1, -1, +2, -2, ... .  Keeping the explicit l <= 4
    // polynomials avoids transcendental functions in the inner grid loop.
    if (l == 1) {
        switch (index) {
            case 0: return  0.4886025119029199f * z;
            case 1: return -0.4886025119029199f * x;
            default:return -0.4886025119029199f * y;
        }
    }

    if (l == 2) {
        switch (index) {
            case 0: return  0.31539156525252005f * (3.0f * z2 - r2);
            case 1: return -1.0925484305920792f * x * z;
            case 2: return -1.0925484305920792f * y * z;
            case 3: return  0.5462742152960396f * (x2 - y2);
            default:return  1.0925484305920792f * x * y;
        }
    }

    if (l == 3) {
        switch (index) {
            case 0: return  0.3731763325901154f * z * (5.0f * z2 - 3.0f * r2);
            case 1: return -0.4570457994644658f * x * (5.0f * z2 - r2);
            case 2: return -0.4570457994644658f * y * (5.0f * z2 - r2);
            case 3: return  1.445305721320277f * z * (x2 - y2);
            case 4: return  2.890611442640554f * x * y * z;
            case 5: return -0.5900435899266435f * x * (x2 - 3.0f * y2);
            default:return -0.5900435899266435f * y * (3.0f * x2 - y2);
        }
    }

    const float z4 = z2 * z2;
    const float r4 = r2 * r2;
    switch (index) {
        case 0:
            return 0.10578554691520431f *
                   (35.0f * z4 - 30.0f * z2 * r2 + 3.0f * r4);
        case 1:
            return -0.6690465435572892f * x * z * (7.0f * z2 - 3.0f * r2);
        case 2:
            return -0.6690465435572892f * y * z * (7.0f * z2 - 3.0f * r2);
        case 3:
            return 0.47308734787878004f * (x2 - y2) * (7.0f * z2 - r2);
        case 4:
            return 0.9461746957575601f * x * y * (7.0f * z2 - r2);
        case 5:
            return -1.7701307697799304f * x * z * (x2 - 3.0f * y2);
        case 6:
            return -1.7701307697799304f * y * z * (3.0f * x2 - y2);
        case 7:
            return 0.6258357354491761f *
                   (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);
        default:
            return 2.5033429417967046f * x * y * (x2 - y2);
    }
}

__device__ __forceinline__ void cartesian_exponents(
    const int l, const int index, int& ax, int& ay, int& az) {
    ax = ay = az = 0;
    if (l == 0) return;

    if (l == 1) {
        const int table[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        ax = table[index][0]; ay = table[index][1]; az = table[index][2];
        return;
    }
    if (l == 2) {
        const int table[6][3] = {
            {2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}
        };
        ax = table[index][0]; ay = table[index][1]; az = table[index][2];
        return;
    }
    if (l == 3) {
        const int table[10][3] = {
            {3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},
            {2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1}
        };
        ax = table[index][0]; ay = table[index][1]; az = table[index][2];
        return;
    }

    const int table[15][3] = {
        {4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},
        {1,3,0},{0,3,1},{1,0,3},{0,1,3},{2,2,0},
        {2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}
    };
    ax = table[index][0]; ay = table[index][1]; az = table[index][2];
}

__device__ __forceinline__ int shell_basis_count_device(const GpuShell& shell) {
    return shell.pure ? (2 * static_cast<int>(shell.l) + 1)
                      : ((static_cast<int>(shell.l) + 1) *
                         (static_cast<int>(shell.l) + 2) / 2);
}

__device__ float evaluate_shell_component(
    const GpuShell& shell,
    const GpuPrimitive* primitives,
    const int component,
    const float dx,
    const float dy,
    const float dz,
    const float r2) {

    if (shell.pure) {
        float radial = 0.0f;
        for (std::uint32_t p = 0; p < shell.primitive_count; ++p) {
            const GpuPrimitive primitive = primitives[shell.primitive_offset + p];
            if (primitive.coefficient == 0.0f) continue;
            const float decay = __expf(-primitive.exponent * r2);
            if (decay == 0.0f) continue;
            const float norm = spherical_primitive_norm(
                primitive.exponent, static_cast<int>(shell.l));
            radial += primitive.coefficient * norm * decay;
        }
        if (radial == 0.0f) return 0.0f;
        return radial * real_solid_harmonic(
            static_cast<int>(shell.l), component, dx, dy, dz);
    }

    int ax, ay, az;
    cartesian_exponents(static_cast<int>(shell.l), component, ax, ay, az);
    float contracted = 0.0f;
    for (std::uint32_t p = 0; p < shell.primitive_count; ++p) {
        const GpuPrimitive primitive = primitives[shell.primitive_offset + p];
        if (primitive.coefficient == 0.0f) continue;
        const float decay = __expf(-primitive.exponent * r2);
        if (decay == 0.0f) continue;
        const float norm = cartesian_primitive_norm(primitive.exponent, ax, ay, az);
        contracted += primitive.coefficient * norm * decay;
    }
    if (contracted == 0.0f) return 0.0f;
    const float monomial = powi(dx, ax) * powi(dy, ay) * powi(dz, az);
    return contracted * monomial;
}

__device__ float orbital_value(
    const GpuShell* shells,
    const std::uint32_t shell_count,
    const GpuPrimitive* primitives,
    const float* coefficients,
    const float x, const float y, const float z) {
    float psi = 0.0f;
    for (std::uint32_t s = 0; s < shell_count; ++s) {
        const GpuShell shell = shells[s];
        const float dx = x - shell.cx;
        const float dy = y - shell.cy;
        const float dz = z - shell.cz;
        const float r2 = dx * dx + dy * dy + dz * dz;

        const int n = shell_basis_count_device(shell);
        for (int c = 0; c < n; ++c) {
            const float coefficient = coefficients[shell.basis_offset + c];
            if (coefficient == 0.0f) continue;
            const float basis = evaluate_shell_component(
                shell, primitives, c, dx, dy, dz, r2);
            psi = fmaf(coefficient, basis, psi);
        }
    }
    return psi;
}

__global__ void orbital_kernel(
    cudaSurfaceObject_t surface,
    const GpuShell* shells,
    const std::uint32_t shell_count,
    const GpuPrimitive* primitives,
    const float* coefficients,
    const GridBox box,
    const int nx,
    const int ny,
    const int nz) {
    const int ix = blockIdx.x * blockDim.x + threadIdx.x;
    const int iy = blockIdx.y * blockDim.y + threadIdx.y;
    const int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= nx || iy >= ny || iz >= nz) return;
    const float tx = nx > 1 ? static_cast<float>(ix) / static_cast<float>(nx - 1) : 0.0f;
    const float ty = ny > 1 ? static_cast<float>(iy) / static_cast<float>(ny - 1) : 0.0f;
    const float tz = nz > 1 ? static_cast<float>(iz) / static_cast<float>(nz - 1) : 0.0f;
    const float x = box.min_x + tx * (box.max_x - box.min_x);
    const float y = box.min_y + ty * (box.max_y - box.min_y);
    const float z = box.min_z + tz * (box.max_z - box.min_z);
    const float psi = orbital_value(shells,shell_count,primitives,coefficients,x,y,z);

    surf3Dwrite(psi, surface,
                static_cast<std::size_t>(ix) * sizeof(float),
                iy, iz);
}

__global__ void orbital_linear_kernel(
    float* output,
    const GpuShell* shells,
    const std::uint32_t shell_count,
    const GpuPrimitive* primitives,
    const float* coefficients,
    const GridBox box,
    const int nx, const int ny, const int nz,
    const int first_z, const int slab_depth) {
    const int ix = blockIdx.x * blockDim.x + threadIdx.x;
    const int iy = blockIdx.y * blockDim.y + threadIdx.y;
    const int local_z = blockIdx.z * blockDim.z + threadIdx.z;
    const int iz = first_z + local_z;
    if (ix >= nx || iy >= ny || local_z >= slab_depth || iz >= nz) return;
    const float tx = nx > 1 ? static_cast<float>(ix) / static_cast<float>(nx - 1) : 0.0f;
    const float ty = ny > 1 ? static_cast<float>(iy) / static_cast<float>(ny - 1) : 0.0f;
    const float tz = nz > 1 ? static_cast<float>(iz) / static_cast<float>(nz - 1) : 0.0f;
    const float x = box.min_x + tx * (box.max_x - box.min_x);
    const float y = box.min_y + ty * (box.max_y - box.min_y);
    const float z = box.min_z + tz * (box.max_z - box.min_z);
    output[static_cast<std::size_t>(ix) + static_cast<std::size_t>(nx) *
        (static_cast<std::size_t>(iy) + static_cast<std::size_t>(ny) * local_z)] =
        orbital_value(shells,shell_count,primitives,coefficients,x,y,z);
}

struct InteropError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

void interop_check(cudaError_t status, const char* what) {
    if (status != cudaSuccess)
        throw InteropError(std::string(what) + ": " + cudaGetErrorString(status));
}

bool gl_uses_device(int device) {
    unsigned int count = 0;
    int devices[16]{};
    const auto status = cudaGLGetDevices(&count,devices,16,cudaGLDeviceListAll);
    if (status != cudaSuccess) {
        (void)cudaGetLastError();
        return false;
    }
    for (unsigned int i=0; i<count; ++i) if (devices[i] == device) return true;
    return false;
}

struct MappedTexture {
    cudaGraphicsResource* resource;
    bool mapped = false;
    explicit MappedTexture(cudaGraphicsResource* value) : resource(value) {
        interop_check(cudaGraphicsMapResources(1,&resource,0),"cudaGraphicsMapResources");
        mapped = true;
    }
    ~MappedTexture() { if (mapped) (void)cudaGraphicsUnmapResources(1,&resource,0); }
    cudaArray_t array() const {
        cudaArray_t result = nullptr;
        interop_check(cudaGraphicsSubResourceGetMappedArray(&result,resource,0,0),
                      "cudaGraphicsSubResourceGetMappedArray");
        return result;
    }
    void unmap() {
        interop_check(cudaGraphicsUnmapResources(1,&resource,0),"cudaGraphicsUnmapResources");
        mapped = false;
    }
};

struct Surface {
    cudaSurfaceObject_t object = 0;
    explicit Surface(cudaArray_t array) {
        cudaResourceDesc desc{};
        desc.resType = cudaResourceTypeArray;
        desc.res.array.array = array;
        const auto status = cudaCreateSurfaceObject(&object,&desc);
        if (status != cudaSuccess) {
            if (object) (void)cudaDestroySurfaceObject(object);
            interop_check(status,"cudaCreateSurfaceObject");
        }
    }
    ~Surface() { if (object) (void)cudaDestroySurfaceObject(object); }
};

struct Events {
    cudaEvent_t start = nullptr;
    cudaEvent_t stop = nullptr;
    Events() {
        const auto status = cudaEventCreate(&start);
        if (status != cudaSuccess) {
            if (start) (void)cudaEventDestroy(start);
            cuda_check(status,"cudaEventCreate start");
        }
        try { cuda_check(cudaEventCreate(&stop),"cudaEventCreate stop"); }
        catch (...) {
            if (stop) (void)cudaEventDestroy(stop);
            (void)cudaEventDestroy(start);
            throw;
        }
    }
    ~Events() {
        if (stop) (void)cudaEventDestroy(stop);
        if (start) (void)cudaEventDestroy(start);
    }
    void begin() { cuda_check(cudaEventRecord(start),"cudaEventRecord start"); }
    double finish() {
        cuda_check(cudaEventRecord(stop),"cudaEventRecord stop");
        cuda_check(cudaEventSynchronize(stop),"cudaEventSynchronize stop");
        float ms = 0.0f;
        cuda_check(cudaEventElapsedTime(&ms,start,stop),"cudaEventElapsedTime");
        return ms;
    }
};

struct DeviceFloatBuffer {
    float* data = nullptr;
    explicit DeviceFloatBuffer(std::size_t count) {
        const auto status = cudaMalloc(&data,count*sizeof(float));
        if (status != cudaSuccess) {
            if (data) (void)cudaFree(data);
            cuda_check(status,"cudaMalloc grid buffer");
        }
    }
    ~DeviceFloatBuffer() { if (data) (void)cudaFree(data); }
};

struct GlUploadState {
    GLint old_binding = 0;
    GLint old_alignment = 4;
    explicit GlUploadState(unsigned int texture) {
        if (!gl::TexSubImage3D) throw std::runtime_error("OpenGL 3D upload function is unavailable");
        glGetIntegerv(GL_TEXTURE_BINDING_3D,&old_binding);
        glGetIntegerv(GL_UNPACK_ALIGNMENT,&old_alignment);
        glBindTexture(GL_TEXTURE_3D,texture);
        glPixelStorei(GL_UNPACK_ALIGNMENT,1);
        while (glGetError() != GL_NO_ERROR) {}
    }
    ~GlUploadState() {
        glPixelStorei(GL_UNPACK_ALIGNMENT,old_alignment);
        glBindTexture(GL_TEXTURE_3D,static_cast<GLuint>(old_binding));
    }
};

} // namespace

struct CudaOrbitalEvaluator::Impl {
    const Wavefunction* wf = nullptr;
    int device = -1;
    GpuShell* d_shells = nullptr;
    GpuPrimitive* d_primitives = nullptr;
    float* d_coefficients = nullptr;
    cudaGraphicsResource* texture_resource = nullptr;
    unsigned int texture = 0;
    std::string device_name_storage = "unknown";
    double last_ms = 0.0;

    explicit Impl(const Wavefunction& wavefunction, int requested_device) : wf(&wavefunction) {
        int device_count = 0;
        cuda_check(cudaGetDeviceCount(&device_count),"cudaGetDeviceCount");
        if (device_count <= 0) throw std::runtime_error("No CUDA device is available");
        if (requested_device < -1 || requested_device >= device_count)
            throw std::out_of_range("Requested CUDA device index is out of range");
        if (requested_device >= 0) device = requested_device;
        else {
            unsigned int gl_count = 0;
            int gl_devices[16]{};
            const auto gl_status = cudaGLGetDevices(&gl_count,gl_devices,16,cudaGLDeviceListAll);
            if (gl_status == cudaSuccess && gl_count > 0 &&
                gl_devices[0] >= 0 && gl_devices[0] < device_count)
                device = gl_devices[0];
            else {
                if (gl_status != cudaSuccess) (void)cudaGetLastError();
                int current = -1;
                if (cudaGetDevice(&current) == cudaSuccess &&
                    current >= 0 && current < device_count) device = current;
                else device = 0;
            }
        }
        cuda_check(cudaSetDevice(device),"cudaSetDevice");

        cudaDeviceProp prop{};
        cuda_check(cudaGetDeviceProperties(&prop, device), "cudaGetDeviceProperties");
        device_name_storage = prop.name;

        std::vector<GpuShell> shells;
        shells.reserve(wavefunction.shells.size());
        for (const Shell& shell : wavefunction.shells) {
            if (shell.angular_momentum > 4 || shell.pure > 1 ||
                shell.primitive_offset > wavefunction.primitives.size() ||
                shell.primitive_count > wavefunction.primitives.size()-shell.primitive_offset ||
                shell.basis_offset > wavefunction.basis_count ||
                shell_basis_count(shell) > wavefunction.basis_count-shell.basis_offset)
                throw std::invalid_argument("Invalid CUDA shell");
            const Atom& atom = wavefunction.atoms.at(shell.atom_index);
            GpuShell gpu{};
            gpu.cx = pack_gpu_scalar(atom.x,"atom x");
            gpu.cy = pack_gpu_scalar(atom.y,"atom y");
            gpu.cz = pack_gpu_scalar(atom.z,"atom z");
            gpu.primitive_offset = shell.primitive_offset;
            gpu.primitive_count = shell.primitive_count;
            gpu.basis_offset = shell.basis_offset;
            gpu.l = shell.angular_momentum;
            gpu.pure = shell.pure;
            shells.push_back(gpu);
        }

        std::vector<GpuPrimitive> primitives;
        primitives.reserve(wavefunction.primitives.size());
        for (const Primitive& p : wavefunction.primitives) {
            if (!std::isfinite(p.exponent) || p.exponent <= 0.0 ||
                !std::isfinite(p.coefficient))
                throw std::invalid_argument("Invalid CUDA Gaussian primitive");
            primitives.push_back({pack_gpu_scalar(p.exponent,"primitive exponent"),
                                  pack_gpu_scalar(p.coefficient,"contraction coefficient")});
        }
        try {
            if (!shells.empty()) {
                cuda_check(cudaMalloc(&d_shells,shells.size()*sizeof(GpuShell)),"cudaMalloc shells");
                cuda_check(cudaMemcpy(d_shells,shells.data(),shells.size()*sizeof(GpuShell),
                                      cudaMemcpyHostToDevice),"cudaMemcpy shells");
            }
            if (!primitives.empty()) {
                cuda_check(cudaMalloc(&d_primitives,primitives.size()*sizeof(GpuPrimitive)),
                           "cudaMalloc primitives");
                cuda_check(cudaMemcpy(d_primitives,primitives.data(),primitives.size()*sizeof(GpuPrimitive),
                                      cudaMemcpyHostToDevice),"cudaMemcpy primitives");
            }
            if (wavefunction.basis_count)
                cuda_check(cudaMalloc(&d_coefficients,
                                      wavefunction.basis_count*sizeof(float)),
                           "cudaMalloc MO coefficients");
        } catch (...) { release(); throw; }
    }

    void bind_device() const { cuda_check(cudaSetDevice(device),"cudaSetDevice"); }
    void release() noexcept {
        if (device < 0 || cudaSetDevice(device) != cudaSuccess) return;
        if (texture_resource) (void)cudaGraphicsUnregisterResource(texture_resource);
        texture_resource = nullptr;
        if (d_coefficients) (void)cudaFree(d_coefficients);
        if (d_primitives) (void)cudaFree(d_primitives);
        if (d_shells) (void)cudaFree(d_shells);
        d_coefficients = nullptr; d_primitives = nullptr; d_shells = nullptr;
    }
    ~Impl() { release(); }
};

CudaOrbitalEvaluator::CudaOrbitalEvaluator(const Wavefunction& wavefunction, int device_index)
    : impl_(std::make_unique<Impl>(wavefunction,device_index)) {}

CudaOrbitalEvaluator::~CudaOrbitalEvaluator() = default;
CudaOrbitalEvaluator::CudaOrbitalEvaluator(CudaOrbitalEvaluator&&) noexcept = default;
CudaOrbitalEvaluator& CudaOrbitalEvaluator::operator=(CudaOrbitalEvaluator&&) noexcept = default;

void CudaOrbitalEvaluator::attach_gl_texture(const unsigned int texture) {
    detach_gl_texture();
    if (!texture) throw std::invalid_argument("CUDA evaluator requires a 3D texture");
    impl_->bind_device();
    impl_->texture = texture;
    if (gl_uses_device(impl_->device)) {
        cudaGraphicsResource* registered = nullptr;
        const auto status = cudaGraphicsGLRegisterImage(
            &registered,texture,GL_TEXTURE_3D,
            cudaGraphicsRegisterFlagsSurfaceLoadStore |
                cudaGraphicsRegisterFlagsWriteDiscard);
        if (status == cudaSuccess) impl_->texture_resource = registered;
        else {
            if (registered) (void)cudaGraphicsUnregisterResource(registered);
            (void)cudaGetLastError();
        }
    }
}

void CudaOrbitalEvaluator::detach_gl_texture() {
    if (impl_ && impl_->texture_resource) {
        impl_->bind_device();
        cuda_check(cudaGraphicsUnregisterResource(impl_->texture_resource),
                   "cudaGraphicsUnregisterResource");
        impl_->texture_resource = nullptr;
    }
    if (impl_) impl_->texture = 0;
}

void CudaOrbitalEvaluator::evaluate(const std::size_t mo_index,
                                    const GridBox& box,
                                    const int nx,
                                    const int ny,
                                    const int nz) {
    if (!impl_->texture) {
        throw std::runtime_error("No OpenGL 3D texture attached to CUDA evaluator");
    }
    if (mo_index >= impl_->wf->orbitals.size()) {
        throw std::out_of_range("MO index out of range");
    }
    if (nx <= 0 || ny <= 0 || nz <= 0 || nx > 512 || ny > 512 || nz > 512) {
        throw std::invalid_argument("CUDA grid dimensions must be between 1 and 512");
    }
    const float bounds[] = {box.min_x,box.min_y,box.min_z,box.max_x,box.max_y,box.max_z};
    for (float value : bounds) if (!std::isfinite(value))
        throw std::invalid_argument("CUDA grid bounds are not finite");
    if (box.min_x > box.max_x || box.min_y > box.max_y || box.min_z > box.max_z ||
        !std::isfinite(box.max_x-box.min_x) ||
        !std::isfinite(box.max_y-box.min_y) ||
        !std::isfinite(box.max_z-box.min_z))
        throw std::invalid_argument("CUDA grid bounds are invalid");
    impl_->bind_device();

    const auto& mo = impl_->wf->orbitals[mo_index];
    if (mo.coefficients.size()!=impl_->wf->basis_count) {
        throw std::runtime_error("MO coefficient dimension does not match GPU basis");
    }
    // Scientific host data remain double. The production float evaluator has
    // an explicit upload representation, validated against its actual texture.
    std::vector<float> gpu_coefficients;
    gpu_coefficients.reserve(mo.coefficients.size());
    for (const double value:mo.coefficients) {
        gpu_coefficients.push_back(pack_gpu_scalar(value,"MO coefficient"));
    }
    if (!gpu_coefficients.empty())
        cuda_check(cudaMemcpy(impl_->d_coefficients,gpu_coefficients.data(),
                              gpu_coefficients.size()*sizeof(float),cudaMemcpyHostToDevice),
                   "cudaMemcpy MO coefficients");

    impl_->last_ms = 0.0;
    const dim3 block(8,8,4);
    if (impl_->texture_resource) {
        try {
            MappedTexture mapped(impl_->texture_resource);
            {
                Surface surface(mapped.array());
                Events events;
                const dim3 grid((nx+block.x-1)/block.x,
                                (ny+block.y-1)/block.y,
                                (nz+block.z-1)/block.z);
                events.begin();
                orbital_kernel<<<grid,block>>>(
                    surface.object,impl_->d_shells,
                    static_cast<std::uint32_t>(impl_->wf->shells.size()),
                    impl_->d_primitives,impl_->d_coefficients,box,nx,ny,nz);
                cuda_check(cudaGetLastError(),"orbital_kernel launch");
                impl_->last_ms = events.finish();
            }
            mapped.unmap();
            return;
        } catch (const InteropError&) {
            // The selected CUDA device can still calculate the grid when the
            // display context cannot share its texture with that device.
            (void)cudaGetLastError();
            cuda_check(cudaGraphicsUnregisterResource(impl_->texture_resource),
                       "cudaGraphicsUnregisterResource after interop failure");
            impl_->texture_resource = nullptr;
        }
    }

    const std::size_t plane = static_cast<std::size_t>(nx)*ny;
    const int slab_limit = static_cast<int>(std::min<std::size_t>(
        nz,std::max<std::size_t>(1,(8u*1024u*1024u)/(plane*sizeof(float)))));
    DeviceFloatBuffer device_output(plane*slab_limit);
    std::vector<float> host_output(plane*slab_limit);
    Events events;
    GlUploadState upload(impl_->texture);
    for (int first_z=0; first_z<nz; first_z+=slab_limit) {
        const int depth = std::min(slab_limit,nz-first_z);
        const dim3 grid((nx+block.x-1)/block.x,
                        (ny+block.y-1)/block.y,
                        (depth+block.z-1)/block.z);
        events.begin();
        orbital_linear_kernel<<<grid,block>>>(
            device_output.data,impl_->d_shells,
            static_cast<std::uint32_t>(impl_->wf->shells.size()),
            impl_->d_primitives,impl_->d_coefficients,box,nx,ny,nz,first_z,depth);
        cuda_check(cudaGetLastError(),"orbital_linear_kernel launch");
        impl_->last_ms += events.finish();
        cuda_check(cudaMemcpy(host_output.data(),device_output.data,
                              plane*depth*sizeof(float),cudaMemcpyDeviceToHost),
                   "cudaMemcpy grid to host");
        gl::TexSubImage3D(GL_TEXTURE_3D,0,0,0,first_z,nx,ny,depth,
                          gl::volume_external_format(),GL_FLOAT,host_output.data());
        if (glGetError() != GL_NO_ERROR)
            throw std::runtime_error("Unable to upload CUDA grid to the display texture");
    }
}

const char* CudaOrbitalEvaluator::device_name() const noexcept {
    return impl_ ? impl_->device_name_storage.c_str() : "unknown";
}

double CudaOrbitalEvaluator::last_kernel_ms() const noexcept {
    return impl_ ? impl_->last_ms : 0.0;
}

} // namespace cov
