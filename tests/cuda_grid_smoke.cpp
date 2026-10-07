// Runtime comparison of the production CUDA evaluator's OpenGL texture with
// the independent CPU grid evaluator. A missing CUDA device is a failure.
#include "cov/cuda_orbital.hpp"
#include "cov/gl_api.hpp"
#include "cov/orbital_grid.hpp"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr GLenum texture_3d = 0x806F;
constexpr GLenum red = 0x1903;
constexpr GLenum r32f = 0x822E;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

struct GlContext {
    GLFWwindow* window = nullptr;
    GlContext() {
        require(glfwInit() == GLFW_TRUE, "GLFW initialization failed");
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        window = glfwCreateWindow(16, 16, "COV CUDA grid smoke", nullptr, nullptr);
        if (!window) {
            glfwTerminate();
            throw std::runtime_error("Hidden OpenGL context creation failed");
        }
        glfwMakeContextCurrent(window);
        if (!cov::gl::load()) {
            glfwDestroyWindow(window);
            window = nullptr;
            glfwTerminate();
            throw std::runtime_error("OpenGL entry points unavailable");
        }
    }
    ~GlContext() {
        if (window) glfwDestroyWindow(window);
        glfwTerminate();
    }
    GlContext(const GlContext&) = delete;
    GlContext& operator=(const GlContext&) = delete;
};

struct Texture {
    GLuint id = 0;
    Texture(int nx, int ny, int nz) {
        glGenTextures(1, &id);
        require(id != 0, "OpenGL texture creation failed");
        glBindTexture(texture_3d, id);
        glTexParameteri(texture_3d, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(texture_3d, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        cov::gl::TexImage3D(texture_3d, 0, r32f, nx, ny, nz, 0, red, GL_FLOAT, nullptr);
        const GLenum error = glGetError();
        glBindTexture(texture_3d, 0);
        if (error != GL_NO_ERROR) {
            glDeleteTextures(1, &id);
            id = 0;
            throw std::runtime_error("OpenGL 3D texture allocation failed");
        }
    }
    ~Texture() { if (id) glDeleteTextures(1, &id); }
    Texture(const Texture&) = delete;
    Texture& operator=(const Texture&) = delete;
};

struct Attachment {
    cov::CudaOrbitalEvaluator& evaluator;
    Attachment(cov::CudaOrbitalEvaluator& evaluator_, GLuint texture)
        : evaluator(evaluator_) { evaluator.attach_gl_texture(texture); }
    ~Attachment() {
        try { evaluator.detach_gl_texture(); } catch (...) {}
    }
    Attachment(const Attachment&) = delete;
    Attachment& operator=(const Attachment&) = delete;
};

struct Case {
    int degree;
    bool pure;
    int component;
    std::size_t basis_index;
};

struct Fixture {
    cov::Wavefunction wavefunction;
    std::vector<Case> cases;
    std::size_t cartesian_g_mo = 0;
    std::size_t pure_g_mo = 0;
    std::size_t zero_mo = 0;
};

Fixture make_fixture() {
    Fixture fixture;
    cov::Atom atom;
    atom.x = 0.2;
    atom.y = -0.3;
    atom.z = 0.1;
    fixture.wavefunction.atoms.push_back(atom);
    for (int degree = 0; degree <= 4; ++degree) {
        for (bool pure : {false, true}) {
            cov::Shell shell;
            shell.primitive_offset = static_cast<std::uint32_t>(fixture.wavefunction.primitives.size());
            shell.primitive_count = 1;
            shell.basis_offset = fixture.wavefunction.basis_count;
            shell.angular_momentum = static_cast<std::uint8_t>(degree);
            shell.pure = pure ? 1 : 0;
            fixture.wavefunction.primitives.push_back({1.25, -0.75});
            fixture.wavefunction.shells.push_back(shell);
            const auto count = cov::shell_basis_count(shell);
            for (std::uint32_t component = 0; component < count; ++component) {
                fixture.cases.push_back({degree, pure, static_cast<int>(component),
                                         static_cast<std::size_t>(shell.basis_offset + component)});
            }
            fixture.wavefunction.basis_count += count;
        }
    }
    for (std::size_t i = 0; i < fixture.cases.size(); ++i) {
        const auto& item = fixture.cases[i];
        cov::MolecularOrbital orbital;
        orbital.coefficients.assign(fixture.wavefunction.basis_count, 0.0);
        orbital.coefficients[item.basis_index] = 0.6;
        fixture.wavefunction.orbitals.push_back(std::move(orbital));
        if (item.degree == 4 && item.component == 0) {
            if (item.pure) fixture.pure_g_mo = i;
            else fixture.cartesian_g_mo = i;
        }
    }
    fixture.zero_mo = fixture.wavefunction.orbitals.size();
    cov::MolecularOrbital zero;
    zero.coefficients.assign(fixture.wavefunction.basis_count, 0.0);
    fixture.wavefunction.orbitals.push_back(std::move(zero));
    return fixture;
}

void compare(cov::CudaOrbitalEvaluator& evaluator, const Fixture& fixture,
             std::size_t mo, const cov::GridBox& box, int nx, int ny, int nz,
             const std::string& label, bool expect_zero = false) {
    const auto request = cov::make_grid_request(box, nx, ny, nz, 0,
                                                cov::checked_grid_size(nx, ny, nz));
    const auto terms = cov::pack_orbital_terms(fixture.wavefunction, mo);
    std::vector<float> reference(request.point_count);
    cov::evaluate_cpu_grid(terms, request, reference);

    Texture texture(nx, ny, nz);
    std::vector<float> actual(request.point_count, std::numeric_limits<float>::quiet_NaN());
    {
        Attachment attachment(evaluator, texture.id);
        evaluator.evaluate(mo, box, nx, ny, nz);
        glBindTexture(texture_3d, texture.id);
        glGetTexImage(texture_3d, 0, red, GL_FLOAT, actual.data());
        const GLenum error = glGetError();
        glBindTexture(texture_3d, 0);
        require(error == GL_NO_ERROR, label + ": OpenGL texture readback failed");
        evaluator.detach_gl_texture();
    }

    float peak = 0.0f;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        require(std::isfinite(reference[i]) && std::isfinite(actual[i]),
                label + ": nonfinite grid value at " + std::to_string(i));
        if (expect_zero) require(reference[i] == 0.0f, label + ": CPU tail is not zero");
        peak = std::max(peak, std::abs(reference[i]));
        const float tolerance = 5.0e-5f * std::max(1.0f, std::abs(reference[i]));
        if (std::abs(actual[i] - reference[i]) > tolerance) {
            throw std::runtime_error(label + ": point " + std::to_string(i) +
                                     " CUDA=" + std::to_string(actual[i]) +
                                     " CPU=" + std::to_string(reference[i]));
        }
    }
    if (!expect_zero) require(peak > 1.0e-4f, label + ": reference has insufficient signal");
}

} // namespace

int main() {
    try {
        GlContext gl;
        const auto fixture = make_fixture();
        cov::CudaOrbitalEvaluator evaluator(fixture.wavefunction, -1);
        const char* device = evaluator.device_name();
        require(device && *device, "CUDA evaluator did not identify a device");
        std::cout << "CUDA device: " << device << '\n';

        const cov::GridBox normal_box{-0.8f, -0.7f, -0.6f, 0.9f, 0.8f, 0.7f};
        for (std::size_t mo = 0; mo < fixture.cases.size(); ++mo) {
            const auto& item = fixture.cases[mo];
            compare(evaluator, fixture, mo, normal_box, 4, 3, 2,
                    "l=" + std::to_string(item.degree) +
                    " pure=" + std::to_string(item.pure) +
                    " component=" + std::to_string(item.component));
        }

        const cov::GridBox tail_box{1.0e10f, 0.0f, 0.0f,
                                    1.0e10f, 0.0f, 0.0f};
        compare(evaluator, fixture, fixture.cartesian_g_mo, tail_box, 1, 1, 1,
                "Cartesian g radial underflow", true);
        compare(evaluator, fixture, fixture.pure_g_mo, tail_box, 1, 1, 1,
                "pure g radial underflow", true);
        compare(evaluator, fixture, fixture.zero_mo, tail_box, 1, 1, 1,
                "zero MO coefficients", true);
        std::cout << "PASS: production CUDA grid smoke\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "FAIL: " << ex.what() << '\n';
        return 1;
    }
}
