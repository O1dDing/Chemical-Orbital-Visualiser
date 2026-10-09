#include "cov/orbital_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

class Module {
public:
    explicit Module(const std::filesystem::path& path) {
        try {
            open(path);
        } catch (...) {
            close();
            throw;
        }
    }

    Module(const Module&) = delete;
    Module& operator=(const Module&) = delete;

    ~Module() { close(); }

    void set_terms(std::span<const CovGaussianTerm> terms) const {
        require(terms.size() <= UINT32_MAX, "Probe has too many terms");
        char error[1024]{};
        const int result = api_->set_terms(context_, terms.data(), static_cast<uint32_t>(terms.size()),
                                          error, sizeof(error));
        require(result != 0, std::string("set_terms failed: ") + error);
    }

    void evaluate(const CovGridRequest& request, std::span<float> values) const {
        require(values.size() >= request.point_count, "Probe output span is too small");
        char error[1024]{};
        const int result = api_->evaluate(context_, &request, values.data(), error, sizeof(error));
        require(result != 0, std::string("evaluate failed: ") + error);
    }

private:
    void open(const std::filesystem::path& path) {
        require(path.is_absolute() && std::filesystem::is_regular_file(path),
                "Pass an existing absolute compute module path");
#ifdef _WIN32
        library_ = LoadLibraryExW(path.c_str(), nullptr,
                                  LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
        if (!library_) throw std::runtime_error("LoadLibraryExW failed: " + std::to_string(GetLastError()));
        auto symbol = GetProcAddress(library_, "cov_get_compute_api");
#else
        library_ = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!library_) throw std::runtime_error(std::string("dlopen failed: ") + dlerror());
        auto symbol = dlsym(library_, "cov_get_compute_api");
#endif
        require(symbol != nullptr, "Compute module has no cov_get_compute_api symbol");
        const auto get_api = reinterpret_cast<CovGetComputeApi>(symbol);
        api_ = get_api(COV_COMPUTE_ABI_VERSION);
        require(api_ && api_->abi_version == COV_COMPUTE_ABI_VERSION &&
                    api_->backend_name && api_->create && api_->destroy &&
                    api_->set_terms && api_->evaluate && api_->device_name,
                "Compute module ABI is missing or incompatible");
        char error[1024]{};
        context_ = api_->create(-1, error, sizeof(error));
        require(context_ != nullptr, std::string("No usable compute device: ") + error);
        const char* name = api_->device_name(context_);
        require(name && *name, "Compute module did not identify its device");
        std::cout << "backend=" << api_->backend_name << " device=" << name << '\n';
    }

    void close() noexcept {
        if (context_) api_->destroy(context_);
        context_ = nullptr;
#ifdef _WIN32
        if (library_) FreeLibrary(library_);
#else
        if (library_) dlclose(library_);
#endif
        library_ = nullptr;
    }
#ifdef _WIN32
    HMODULE library_ = nullptr;
#else
    void* library_ = nullptr;
#endif
    const CovComputeApi* api_ = nullptr;
    void* context_ = nullptr;
};

void compare(Module& module, const std::string& label,
             std::span<const CovGaussianTerm> terms, const CovGridRequest& request,
             bool exact_zero = false) {
    module.set_terms(terms);
    std::vector<float> reference(request.point_count);
    std::vector<float> actual(request.point_count, std::numeric_limits<float>::quiet_NaN());
    cov::evaluate_cpu_grid(terms, request, reference);
    module.evaluate(request, actual);
    for (size_t i = 0; i < actual.size(); ++i) {
        require(std::isfinite(reference[i]) && std::isfinite(actual[i]),
                label + ": nonfinite result at point " + std::to_string(i));
        if (exact_zero) require(reference[i] == 0.0f, label + ": analytic zero disagrees with CPU");
        const float tolerance = 5.0e-5f * std::max(1.0f, std::abs(reference[i]));
        if (std::abs(actual[i] - reference[i]) > tolerance) {
            throw std::runtime_error(label + ": point " + std::to_string(i) +
                                     " device=" + std::to_string(actual[i]) +
                                     " CPU=" + std::to_string(reference[i]));
        }
    }
}

cov::Wavefunction one_shell(int degree, bool pure, int component) {
    cov::Wavefunction wf;
    cov::Atom atom;
    atom.x = 0.2;
    atom.y = -0.3;
    atom.z = 0.1;
    wf.atoms.push_back(atom);
    wf.primitives.push_back({1.25, -0.75});
    cov::Shell shell;
    shell.primitive_count = 1;
    shell.angular_momentum = static_cast<std::uint8_t>(degree);
    shell.pure = pure ? 1 : 0;
    wf.shells.push_back(shell);
    wf.basis_count = cov::shell_basis_count(shell);
    cov::MolecularOrbital mo;
    mo.coefficients.assign(wf.basis_count, 0.0);
    mo.coefficients[component] = 0.6;
    wf.orbitals.push_back(std::move(mo));
    return wf;
}

void packed_components(Module& module) {
    const cov::GridBox point{0.63f, -0.67f, 0.69f, 0.63f, -0.67f, 0.69f};
    const auto request = cov::make_grid_request(point, 1, 1, 1, 0, 1);
    for (int degree = 0; degree <= 4; ++degree) {
        for (bool pure : {false, true}) {
            const int components = pure ? 2 * degree + 1 : (degree + 1) * (degree + 2) / 2;
            for (int component = 0; component < components; ++component) {
                const auto terms = cov::pack_orbital_terms(one_shell(degree, pure, component), 0);
                compare(module, "packed l=" + std::to_string(degree) +
                                " pure=" + std::to_string(pure) +
                                " component=" + std::to_string(component), terms, request);
            }
        }
    }
}

std::vector<CovGaussianTerm> contracted_terms() {
    cov::Wavefunction wf;
    cov::Atom first, second;
    first.x = -0.4; first.y = 0.1; first.z = 0.2;
    second.x = 0.7; second.y = -0.2; second.z = -0.1;
    wf.atoms = {first, second};
    wf.primitives = {{0.8, 0.7}, {2.2, -0.2}, {0.6, 0.9}, {1.7, 0.15}};
    cov::Shell s, p;
    s.atom_index = 0; s.primitive_offset = 0; s.primitive_count = 2;
    s.basis_offset = 0; s.angular_momentum = 0;
    p.atom_index = 1; p.primitive_offset = 2; p.primitive_count = 2;
    p.basis_offset = 1; p.angular_momentum = 1;
    wf.shells = {s, p};
    wf.basis_count = 4;
    cov::MolecularOrbital mo;
    mo.coefficients = {0.55, -0.35, 0.2, 0.8};
    wf.orbitals.push_back(std::move(mo));
    return cov::pack_orbital_terms(wf, 0);
}

void grids_and_chunks(Module& module) {
    const auto terms = contracted_terms();
    const cov::GridBox box{-0.8f, -0.5f, -0.6f, 1.0f, 0.8f, 0.7f};
    constexpr std::array<std::array<int, 3>, 4> shapes{{{5, 3, 2}, {1, 4, 3},
                                                         {4, 1, 3}, {4, 3, 1}}};
    for (const auto& dims : shapes) {
        const auto count = cov::checked_grid_size(dims[0], dims[1], dims[2]);
        for (size_t first = 0; first < count; first += 7) {
            const auto n = std::min<size_t>(7, count - first);
            compare(module, "contracted chunk " + std::to_string(dims[0]) + "x" +
                            std::to_string(dims[1]) + "x" + std::to_string(dims[2]) +
                            " at " + std::to_string(first),
                    terms, cov::make_grid_request(box, dims[0], dims[1], dims[2], first, n));
        }
    }
}

void cancellation_and_tail(Module& module) {
    const cov::GridBox origin{0.25f, 0.0f, 0.0f, 0.25f, 0.0f, 0.0f};
    const auto near_request = cov::make_grid_request(origin, 1, 1, 1, 0, 1);
    const std::vector<CovGaussianTerm> cancelling{
        {0, 0, 0, 1, 1, 0, 0, 0},
        {0, 0, 0, 1, -1, 0, 0, 0},
        {0, 0, 0, 1, 0.125f, 0, 0, 0}
    };
    compare(module, "cancellation", cancelling, near_request);

    const cov::GridBox far_box{1.0e10f, 0, 0, 1.0e10f, 0, 0};
    const auto far_request = cov::make_grid_request(far_box, 1, 1, 1, 0, 1);
    const std::vector<CovGaussianTerm> underflow{
        {0, 0, 0, 1, 1, 4, 0, 0}
    };
    compare(module, "g radial underflow at x=1e10", underflow, far_request, true);

    const std::vector<CovGaussianTerm> zero_coefficient{
        {0, 0, 0, std::numeric_limits<float>::denorm_min(), 0, 4, 0, 0}
    };
    compare(module, "zero coefficient with overflowing power", zero_coefficient, far_request, true);
}

} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "Usage: cov_compute_module_probe ABSOLUTE_MODULE_PATH");
        Module module{std::filesystem::path(argv[1])};
        packed_components(module);
        grids_and_chunks(module);
        cancellation_and_tail(module);
        std::cout << "PASS: native compute module numerical conformance\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "FAIL: " << ex.what() << '\n';
        return 1;
    }
}
