#include "cov/orbital_evaluator.hpp"
#include "cov/volume_renderer.hpp"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
constexpr GLenum texture_3d = 0x806f;
constexpr double pi = 3.14159265358979323846;

void require(bool condition, const char* why) {
    if (!condition) throw std::runtime_error(why);
}

struct Window {
    GLFWwindow* handle = nullptr;
    Window() {
        if (!glfwInit()) throw std::runtime_error("GLFW initialization failed");
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
        handle = glfwCreateWindow(64, 64, "COV compute smoke", nullptr, nullptr);
        if (!handle) { glfwTerminate(); throw std::runtime_error("OpenGL context creation failed"); }
        glfwMakeContextCurrent(handle);
    }
    ~Window() {
        if (handle) glfwDestroyWindow(handle);
        glfwTerminate();
    }
};

cov::Wavefunction sample() {
    cov::Wavefunction wf;
    cov::Atom atom;
    atom.symbol = "H";
    atom.atomic_number = 1;
    atom.x = 0.25;
    atom.y = -0.5;
    atom.z = 0.125;
    wf.atoms.push_back(atom);
    wf.primitives.push_back({1.25, -0.75});
    cov::Shell shell;
    shell.primitive_count = 1;
    wf.shells.push_back(shell);
    wf.basis_count = 1;
    cov::MolecularOrbital plus, minus;
    plus.coefficients = {0.6};
    minus.coefficients = {-0.6};
    wf.orbitals.push_back(plus);
    wf.orbitals.push_back(minus);
    return wf;
}

std::vector<float> read_volume(const cov::VolumeRenderer& renderer) {
    const auto count = std::size_t(renderer.nx()) * renderer.ny() * renderer.nz();
    std::vector<float> values(count);
    glBindTexture(texture_3d, renderer.volume_texture());
    glGetTexImage(texture_3d, 0, GL_RED, GL_FLOAT, values.data());
    const GLenum error = glGetError();
    glBindTexture(texture_3d, 0);
    require(error == GL_NO_ERROR, "OpenGL texture readback failed");
    return values;
}

void expected_volume(const std::vector<float>& values, const cov::GridBox& box,
                     int nx, int ny, int nz, double orbital_sign) {
    require(values.size() == std::size_t(nx) * ny * nz, "Unexpected grid size");
    const double norm = std::pow(2.0 * 1.25 / pi, 0.75);
    for (int z = 0; z < nz; ++z) for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
        const double px = box.min_x + (box.max_x - box.min_x) * (nx == 1 ? 0.0 : double(x)/(nx-1));
        const double py = box.min_y + (box.max_y - box.min_y) * (ny == 1 ? 0.0 : double(y)/(ny-1));
        const double pz = box.min_z + (box.max_z - box.min_z) * (nz == 1 ? 0.0 : double(z)/(nz-1));
        const double dx = px - 0.25, dy = py + 0.5, dz = pz - 0.125;
        const double reference = orbital_sign * 0.6 * -0.75 * norm *
            std::exp(-1.25 * (dx*dx + dy*dy + dz*dz));
        const float actual = values[std::size_t(x) + std::size_t(nx) * (y + ny*z)];
        const double tolerance = 3.0e-5 * std::max(1.0, std::abs(reference));
        if (!std::isfinite(actual) || std::abs(actual - reference) > tolerance)
            throw std::runtime_error("OpenGL volume differs from analytic s orbital");
    }
}

void finish(cov::OrbitalEvaluator& evaluator) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(15);
    while (!evaluator.poll()) {
        require(std::chrono::steady_clock::now() < deadline, "Orbital evaluation timed out");
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    require(evaluator.ready(), "Finished orbital is not ready");
}

void cpu_texture_and_cancel(const cov::Wavefunction& wf, cov::VolumeRenderer& renderer) {
    cov::OrbitalEvaluator evaluator(wf, {"cpu", -1});
    evaluator.attach_gl_texture(renderer.volume_texture());
    const cov::GridBox box{-1.0f, -0.5f, -0.75f, 1.0f, 0.5f, 0.75f};
    renderer.resize_volume(3, 1, 2); // non-cubic; y is a singleton dimension
    evaluator.evaluate(0, box, 3, 1, 2);
    require(evaluator.ready(), "Synchronous CPU evaluation is not ready");
    const auto initial = read_volume(renderer);
    expected_volume(initial, box, 3, 1, 2, 1.0);

    // Cancellation after begin must never publish its pending result.
    evaluator.begin_evaluate(1, box, 3, 1, 2);
    evaluator.cancel();
    require(!evaluator.ready() && !evaluator.poll(), "Cancelled result was published");
    require(read_volume(renderer) == initial, "Cancellation changed the displayed texture");

    // Starting a new MO joins any old worker before the replacement result is uploaded.
    evaluator.begin_evaluate(0, box, 3, 1, 2);
    evaluator.begin_evaluate(1, box, 3, 1, 2);
    finish(evaluator);
    expected_volume(read_volume(renderer), box, 3, 1, 2, -1.0);
}

void auto_and_explicit(const cov::Wavefunction& wf, cov::VolumeRenderer& renderer,
                       const std::string& mode) {
    const cov::GridBox box{-1.0f, -0.5f, -0.75f, 1.0f, 0.5f, 0.75f};
    renderer.resize_volume(3, 1, 2);
    if (mode == "missing") {
        cov::OrbitalEvaluator automatic(wf, {"auto", -1});
        require(std::string(automatic.device_name()) == "CPU", "Missing modules did not select CPU");
        automatic.attach_gl_texture(renderer.volume_texture());
        automatic.evaluate(0, box, 3, 1, 2);
        expected_volume(read_volume(renderer), box, 3, 1, 2, 1.0);
        bool failed = false;
        try { cov::OrbitalEvaluator explicit_hip(wf, {"hip", -1}); }
        catch (const std::exception& ex) {
            failed = std::string(ex.what()).find("hip compute module is unavailable") != std::string::npos;
        }
        require(failed, "Missing explicit HIP backend did not report module absence");
        return;
    }
    require(mode == "failing", "Unknown smoke test mode");
    cov::OrbitalEvaluator automatic(wf, {"auto", -1});
    require(std::string(automatic.device_name()).find("fake HIP GPU") != std::string::npos,
            "Auto did not select the fake HIP module");
    automatic.attach_gl_texture(renderer.volume_texture());
    automatic.evaluate(0, box, 3, 1, 2);
    require(std::string(automatic.device_name()) == "CPU", "Failed HIP did not fall back to CPU");
    expected_volume(read_volume(renderer), box, 3, 1, 2, 1.0);

    cov::OrbitalEvaluator explicit_hip(wf, {"hip", -1});
    explicit_hip.attach_gl_texture(renderer.volume_texture());
    bool failed = false;
    try { explicit_hip.evaluate(0, box, 3, 1, 2); }
    catch (const std::exception& ex) {
        failed = std::string(ex.what()).find("fake HIP evaluate failure") != std::string::npos;
    }
    require(failed, "Explicit HIP error was lost or silently fell back");
}

void native_texture(const cov::Wavefunction& wf, cov::VolumeRenderer& renderer,
                    const std::string& backend) {
    const auto* original_version = glGetString(GL_VERSION);
    require(original_version != nullptr, "No current OpenGL context before module creation");
    const std::string version(reinterpret_cast<const char*>(original_version));
    const cov::GridBox box{-1.0f, -0.5f, -0.75f, 1.0f, 0.5f, 0.75f};
    {
        cov::OrbitalEvaluator evaluator(wf, {backend, -1});
        const auto* current_version = glGetString(GL_VERSION);
        require(current_version && version == reinterpret_cast<const char*>(current_version),
                "Compute module replaced or cleared the viewer OpenGL context");
        require(std::string(evaluator.device_name()).starts_with(backend),
                "Explicit native test did not use the requested backend");
        std::cout << "device=" << evaluator.device_name() << '\n';
        renderer.resize_volume(3, 1, 2);
        evaluator.attach_gl_texture(renderer.volume_texture());
        evaluator.evaluate(0, box, 3, 1, 2);
        expected_volume(read_volume(renderer), box, 3, 1, 2, 1.0);
        evaluator.begin_evaluate(1, box, 3, 1, 2);
        finish(evaluator);
        expected_volume(read_volume(renderer), box, 3, 1, 2, -1.0);
    }
    require(glGetString(GL_VERSION) != nullptr, "Module teardown cleared the viewer OpenGL context");
    expected_volume(read_volume(renderer), box, 3, 1, 2, -1.0);
}
}

int main(int argc, char** argv) {
    try {
        require(argc == 2, "Expected missing, failing or native backend (webgpu/hip/sycl/metal) mode");
        Window window;
        cov::VolumeRenderer renderer;
        const auto wf = sample();
        cpu_texture_and_cancel(wf, renderer);
        const std::string mode = argv[1];
        if (mode == "webgpu" || mode == "hip" || mode == "sycl" || mode == "metal")
            native_texture(wf, renderer, mode);
        else
            auto_and_explicit(wf, renderer, mode);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "viewer_compute_smoke: " << ex.what() << '\n';
        return 1;
    }
}
