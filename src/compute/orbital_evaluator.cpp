#include "cov/orbital_evaluator.hpp"
#include "cov/orbital_grid.hpp"
#include "cov/gl_api.hpp"
#include "cov/threading.hpp"
#ifdef COV_ENABLE_CUDA
#include "cov/cuda_orbital.hpp"
#endif
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>
#ifndef _WIN32
#include <dlfcn.h>
#include <unistd.h>
#endif
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#ifdef __FreeBSD__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

namespace cov {
namespace {
constexpr GLenum texture_3d = 0x806F;

std::filesystem::path module_directory() {
    if (const char* override_dir = std::getenv("COV_COMPUTE_MODULE_DIR")) {
        auto path = std::filesystem::u8path(override_dir);
        if (!path.is_absolute()) throw std::runtime_error("COV_COMPUTE_MODULE_DIR must be an absolute path");
        return path;
    }
#ifdef _WIN32
    std::wstring path(32768, L'\0');
    const auto size = GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
    if (!size || size == path.size()) throw std::runtime_error("Unable to locate the executable directory");
    path.resize(size);
#elif defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::string path(size, '\0');
    if (_NSGetExecutablePath(path.data(), &size) != 0) throw std::runtime_error("Unable to locate the executable directory");
    path.resize(path.find('\0'));
#elif defined(__FreeBSD__)
    int mib[] = {CTL_KERN, KERN_PROC, KERN_PROC_PATHNAME, -1};
    std::string path(4096, '\0');
    size_t size = path.size();
    if (sysctl(mib, 4, path.data(), &size, nullptr, 0) != 0) throw std::runtime_error("Unable to locate the executable directory");
    path.resize(path.find('\0'));
#else
    std::string path(4096, '\0');
    const auto size = readlink("/proc/self/exe", path.data(), path.size());
    if (size <= 0 || static_cast<size_t>(size) == path.size()) throw std::runtime_error("Unable to locate the executable directory");
    path.resize(static_cast<size_t>(size));
#endif
    return std::filesystem::canonical(path).parent_path();
}

struct Module {
#ifdef _WIN32
    HMODULE library = nullptr;
    DLL_DIRECTORY_COOKIE runtime_directory = nullptr;
#else
    void* library = nullptr;
#endif
    const CovComputeApi* api = nullptr;
    void* context = nullptr;
    ~Module() {
        if (context) api->destroy(context);
#ifdef _WIN32
        if (library) FreeLibrary(library);
        if (runtime_directory) RemoveDllDirectory(runtime_directory);
#else
        if (library) dlclose(library);
#endif
    }
    explicit Module(const std::string& name, int device) {
#ifdef _WIN32
        const auto path = module_directory() / ("cov_compute_" + name + ".dll");
        if (const char* runtime = std::getenv("COV_COMPUTE_RUNTIME_DIR")) {
            const auto runtime_path = std::filesystem::u8path(runtime);
            if (!runtime_path.is_absolute()) throw std::runtime_error("COV_COMPUTE_RUNTIME_DIR must be an absolute path");
            runtime_directory = AddDllDirectory(runtime_path.c_str());
            if (!runtime_directory) throw std::runtime_error("Unable to add the compute runtime directory");
        }
        library = LoadLibraryExW(path.c_str(), nullptr, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
        auto get = library ? reinterpret_cast<CovGetComputeApi>(GetProcAddress(library, "cov_get_compute_api")) : nullptr;
#else
#ifdef __APPLE__
        constexpr auto extension = ".dylib";
#else
        constexpr auto extension = ".so";
#endif
        const auto path = module_directory() / ("libcov_compute_" + name + extension);
        library = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        auto get = library ? reinterpret_cast<CovGetComputeApi>(dlsym(library, "cov_get_compute_api")) : nullptr;
#endif
        // A throwing constructor cannot run its destructor; release the library
        // explicitly if ABI negotiation or device creation fails.
        try {
            if (!get) throw std::runtime_error(name + " compute module is unavailable");
            api = get(COV_COMPUTE_ABI_VERSION);
            if (!api || api->abi_version != COV_COMPUTE_ABI_VERSION || !api->create || !api->destroy ||
                !api->set_terms || !api->evaluate || !api->device_name)
                throw std::runtime_error(name + " compute module has an incompatible interface");
            char error[1024]{};
            context = api->create(device, error, sizeof(error));
            if (!context) throw std::runtime_error(name + ": " + error);
        } catch (...) {
#ifdef _WIN32
            if (library) FreeLibrary(library);
            if (runtime_directory) RemoveDllDirectory(runtime_directory);
            runtime_directory = nullptr;
#else
            if (library) dlclose(library);
#endif
            library = nullptr;
            throw;
        }
    }
};
}

struct OrbitalEvaluator::Impl {
    const Wavefunction& wavefunction;
    ComputeOptions options;
    std::unique_ptr<Module> module;
#ifdef COV_ENABLE_CUDA
    std::unique_ptr<CudaOrbitalEvaluator> cuda;
#endif
    std::string device = "CPU";
    unsigned texture = 0;
    int nx = 0, ny = 0, nz = 0;
    bool valid = false;
    bool running = false;
    double elapsed = 0;
    std::vector<float> result;
    std::exception_ptr failure;
    std::atomic<bool> done{false};
    cov::jthread worker; // destroyed before the module and wavefunction reference

    Impl(const Wavefunction& wf, ComputeOptions o) : wavefunction(wf), options(std::move(o)) {
        const std::vector<std::string> supported{"auto", "cpu", "cuda", "hip", "sycl", "metal", "webgpu"};
        if (std::find(supported.begin(), supported.end(), options.backend) == supported.end())
            throw std::invalid_argument("Unknown compute backend: " + options.backend);
        if (options.device_index < -1) throw std::invalid_argument("Compute device index must be nonnegative");
        if (options.backend == "cpu") {
            if (options.device_index != -1) throw std::invalid_argument("CPU does not accept a GPU device index");
            return;
        }
        if (options.backend != "auto") { select(options.backend); return; }
        const auto vendor_ptr = glGetString(GL_VENDOR);
        const std::string vendor = vendor_ptr ? reinterpret_cast<const char*>(vendor_ptr) : "";
        std::vector<std::string> candidates;
#ifdef __APPLE__
        candidates = {"metal"};
#else
        if (vendor.find("NVIDIA") != std::string::npos) candidates.push_back("cuda");
        if (vendor.find("AMD") != std::string::npos || vendor.find("ATI") != std::string::npos) candidates.push_back("hip");
        if (vendor.find("Intel") != std::string::npos) candidates.push_back("sycl");
        for (const auto* candidate : {"cuda", "hip", "sycl"})
            if (std::find(candidates.begin(), candidates.end(), candidate) == candidates.end()) candidates.emplace_back(candidate);
#endif
        for (const auto& candidate : candidates) {
            try { select(candidate); return; } catch (const std::exception&) {}
        }
        if (options.device_index >= 0) throw std::runtime_error("No native compute backend accepts the requested device index");
    }
    void select(const std::string& backend) {
        if (backend == "cuda") {
#ifdef COV_ENABLE_CUDA
            cuda = std::make_unique<CudaOrbitalEvaluator>(wavefunction, options.device_index);
            device = std::string("CUDA · ") + cuda->device_name();
            return;
#else
            throw std::runtime_error("This build does not include CUDA");
#endif
        }
        auto selected = std::make_unique<Module>(backend, options.device_index);
        device = backend + " · " + selected->api->device_name(selected->context);
        module = std::move(selected);
    }
    void cancel() {
        if (worker.joinable()) { worker.request_stop(); worker.join(); }
        running = false;
        done.store(false);
        failure = nullptr;
        result.clear();
    }
    bool allow_fallback() const { return options.backend == "auto" && options.device_index == -1; }
    void use_cpu() {
        module.reset();
#ifdef COV_ENABLE_CUDA
        cuda.reset();
#endif
        device = "CPU";
    }
    void upload() {
        GLint old_binding = 0, old_alignment = 4;
        glGetIntegerv(0x806A /* GL_TEXTURE_BINDING_3D */, &old_binding);
        glGetIntegerv(GL_UNPACK_ALIGNMENT, &old_alignment);
        glBindTexture(texture_3d, texture);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        while (glGetError() != GL_NO_ERROR) {}
        gl::TexSubImage3D(texture_3d, 0, 0, 0, 0, nx, ny, nz, gl::volume_external_format(), GL_FLOAT, result.data());
        const auto error = glGetError();
        glPixelStorei(GL_UNPACK_ALIGNMENT, old_alignment);
        glBindTexture(texture_3d, static_cast<GLuint>(old_binding));
        if (error != GL_NO_ERROR) throw std::runtime_error("Unable to upload the orbital grid to the display texture");
        valid = true;
    }
};

OrbitalEvaluator::OrbitalEvaluator(const Wavefunction& wf, ComputeOptions options)
    : impl_(std::make_unique<Impl>(wf, std::move(options))) {}
OrbitalEvaluator::~OrbitalEvaluator() { impl_->cancel(); }
void OrbitalEvaluator::cancel() { impl_->cancel(); impl_->valid = false; }
void OrbitalEvaluator::attach_gl_texture(unsigned texture) {
    detach_gl_texture();
    impl_->texture = texture;
#ifdef COV_ENABLE_CUDA
    if (impl_->cuda) {
        try { impl_->cuda->attach_gl_texture(texture); }
        catch (...) { if (!impl_->allow_fallback()) throw; impl_->use_cpu(); }
    }
#endif
}
void OrbitalEvaluator::detach_gl_texture() {
    cancel();
#ifdef COV_ENABLE_CUDA
    if (impl_->cuda) impl_->cuda->detach_gl_texture();
#endif
    impl_->texture = 0;
}
void OrbitalEvaluator::begin_evaluate(std::size_t mo, const GridBox& box, int nx, int ny, int nz) {
    cancel();
    if (!impl_->texture) throw std::runtime_error("No display texture is attached");
    const auto size = checked_grid_size(nx, ny, nz);
    (void)make_grid_request(box, nx, ny, nz, 0, size);
    impl_->nx = nx; impl_->ny = ny; impl_->nz = nz;
#ifdef COV_ENABLE_CUDA
    if (impl_->cuda) {
        try {
            impl_->cuda->evaluate(mo, box, nx, ny, nz);
            impl_->elapsed = impl_->cuda->last_kernel_ms();
            impl_->valid = true;
            return;
        } catch (...) { if (!impl_->allow_fallback()) throw; impl_->use_cpu(); }
    }
#endif
    auto terms = pack_orbital_terms(impl_->wavefunction, mo);
    if (terms.size() > std::numeric_limits<uint32_t>::max()) throw std::runtime_error("Orbital basis is too large");
    impl_->result.resize(size);
    impl_->running = true;
    impl_->worker = cov::jthread([p = impl_.get(), terms = std::move(terms), box, size](cov::stop_token stop) {
        try {
            const auto start = std::chrono::steady_clock::now();
            auto compute = [&] {
                char error[1024]{};
                if (p->module && !p->module->api->set_terms(p->module->context, terms.data(), static_cast<uint32_t>(terms.size()), error, sizeof(error)))
                    throw std::runtime_error(error);
                constexpr std::size_t chunk_size = 262144;
                for (std::size_t first = 0; first < size; first += chunk_size) {
                    if (stop.stop_requested()) return;
                    const auto count = std::min(chunk_size, size - first);
                    const auto request = make_grid_request(box, p->nx, p->ny, p->nz, first, count);
                    auto output = std::span<float>(p->result).subspan(first, count);
                    if (p->module) {
                        if (!p->module->api->evaluate(p->module->context, &request, output.data(), error, sizeof(error)))
                            throw std::runtime_error(error);
                        if (!std::all_of(output.begin(), output.end(), [](float value) { return std::isfinite(value); }))
                            throw std::runtime_error("The compute device returned a non-finite orbital grid");
                    } else evaluate_cpu_grid(terms, request, output, stop);
                }
            };
            try { compute(); }
            catch (...) {
                if (!p->module || !p->allow_fallback() || stop.stop_requested()) throw;
                // Module state is worker-owned until join; the UI sees the new
                // device label only after the completion acquire below.
                p->module.reset();
                compute();
            }
            p->elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        } catch (...) { p->failure = std::current_exception(); }
        p->done.store(true, std::memory_order_release);
    });
}
bool OrbitalEvaluator::poll() {
    if (!impl_->running || !impl_->done.load(std::memory_order_acquire)) return false;
    if (impl_->worker.joinable()) impl_->worker.join();
    impl_->running = false;
    if (!impl_->module) impl_->device = "CPU";
    if (impl_->failure) std::rethrow_exception(impl_->failure);
    impl_->upload();
    impl_->result.clear();
    return true;
}
void OrbitalEvaluator::evaluate(std::size_t mo, const GridBox& box, int nx, int ny, int nz) {
    begin_evaluate(mo, box, nx, ny, nz);
    if (impl_->running) { impl_->worker.join(); poll(); }
}
bool OrbitalEvaluator::busy() const noexcept { return impl_->running; }
bool OrbitalEvaluator::ready() const noexcept { return impl_->valid; }
const char* OrbitalEvaluator::device_name() const noexcept { return impl_->device.c_str(); }
double OrbitalEvaluator::last_kernel_ms() const noexcept { return impl_->running ? 0.0 : impl_->elapsed; }
} // namespace cov
