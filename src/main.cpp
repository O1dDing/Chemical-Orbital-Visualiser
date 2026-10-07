#include "cov/orbital_evaluator.hpp"
#include "cov/numerical_diagnostics.hpp"
#include "cov/pi_topology_evidence.hpp"
#include <sstream>
#include "cov/file_dialog.hpp"
#include "cov/gl_api.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/molden_parser.hpp"
#include "cov/molecule_style.hpp"
#include "cov/orbital_tracking.hpp"
#include "cov/orbital_ui.hpp"
#include "cov/ui.hpp"
#include "cov/volume_renderer.hpp"
#include "cov/validation.hpp"
#include "cov/viewer_layout.hpp"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <shellapi.h>
#endif

namespace {

std::string g_dropped_path;

enum class StatusKind {
    Ready,
    Parsing,
    Loaded,
    GridUpdated,
    Computing,
    Exported,
    Error,
};

void drop_callback(GLFWwindow*, const int count, const char** paths) {
    if (count > 0 && paths && paths[0]) {
        g_dropped_path = paths[0];
    }
}

cov::GridBox make_grid_box(const cov::Wavefunction& wf, const float padding_bohr = 4.0f) {
    if (wf.atoms.empty()) return {};

    float min_x = static_cast<float>(wf.atoms.front().x);
    float min_y = static_cast<float>(wf.atoms.front().y);
    float min_z = static_cast<float>(wf.atoms.front().z);
    float max_x = min_x;
    float max_y = min_y;
    float max_z = min_z;

    for (const auto& atom : wf.atoms) {
        min_x = std::min(min_x, static_cast<float>(atom.x));
        min_y = std::min(min_y, static_cast<float>(atom.y));
        min_z = std::min(min_z, static_cast<float>(atom.z));
        max_x = std::max(max_x, static_cast<float>(atom.x));
        max_y = std::max(max_y, static_cast<float>(atom.y));
        max_z = std::max(max_z, static_cast<float>(atom.z));
    }

    const float cx = 0.5f * (min_x + max_x);
    const float cy = 0.5f * (min_y + max_y);
    const float cz = 0.5f * (min_z + max_z);
    const float extent = std::max({max_x - min_x, max_y - min_y, max_z - min_z});
    const float half = 0.5f * extent + padding_bohr;

    return {cx-half, cy-half, cz-half, cx+half, cy+half, cz+half};
}

std::size_t initial_orbital(const cov::Wavefunction& wf) {
    const auto frontier=cov::find_frontier_orbitals(wf.orbitals);
    return frontier.homo.value_or(0u);
}

const char* status_label(const StatusKind status, const cov::ui::Language language) {
    using cov::ui::Text;
    switch (status) {
        case StatusKind::Computing:
            switch (language) {
                case cov::ui::Language::ChineseSimplified: return "正在计算轨道网格";
                case cov::ui::Language::Japanese: return "軌道グリッドを計算中";
                case cov::ui::Language::French: return "Calcul de la grille orbitale";
                default: return "Computing orbital grid";
            }
        case StatusKind::Parsing: return cov::ui::tr(Text::Parsing, language);
        case StatusKind::Loaded: return cov::ui::tr(Text::Loaded, language);
        case StatusKind::GridUpdated: return cov::ui::tr(Text::GridUpdated, language);
        case StatusKind::Exported: return cov::ui::tr(Text::Exported, language);
        case StatusKind::Error: return cov::ui::tr(Text::Error, language);
        default: return cov::ui::tr(Text::Ready, language);
    }
}

cov::ui::Tone status_tone(const StatusKind status) {
    switch (status) {
        case StatusKind::Loaded:
        case StatusKind::GridUpdated:
        case StatusKind::Exported: return cov::ui::Tone::Success;
        case StatusKind::Error: return cov::ui::Tone::Danger;
        case StatusKind::Parsing: return cov::ui::Tone::Accent;
        default: return cov::ui::Tone::Neutral;
    }
}

std::filesystem::path path_from_utf8(const std::string& value) {
#ifdef _WIN32
    return std::filesystem::u8path(value);
#else
    return std::filesystem::path(value);
#endif
}

std::string path_to_utf8(const std::filesystem::path& path) {
#if defined(__cpp_lib_char8_t)
    const auto value = path.u8string();
    return std::string(reinterpret_cast<const char*>(value.data()), value.size());
#else
    return path.u8string();
#endif
}

void copy_path_to_buffer(const std::filesystem::path& path,
                         std::array<char, 2048>& buffer) {
    const std::string value = path_to_utf8(path);
    std::snprintf(buffer.data(), buffer.size(), "%s", value.c_str());
}

void disabled_wrapped(const char* text) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    ImGui::TextWrapped("%s", text);
    ImGui::PopStyleColor();
}

void metric_row(const char* label, const char* value) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    disabled_wrapped(label);
    ImGui::TableNextColumn();
    ImGui::TextWrapped("%s", value);
}

void push_recent(std::vector<std::filesystem::path>& recent,
                 const std::filesystem::path& path) {
    const auto normalized = path.lexically_normal();
    recent.erase(std::remove_if(recent.begin(), recent.end(), [&](const auto& existing) {
        return existing.lexically_normal() == normalized;
    }), recent.end());
    recent.insert(recent.begin(), normalized);
    if (recent.size() > 8) recent.resize(8);
}

const char* molecule_style_name(const cov::MoleculeStyle style,
                                const cov::ui::Language language) {
    return cov::ui::tr(style == cov::MoleculeStyle::StickDelocalisation
                           ? cov::ui::Text::StickDelocalisation
                           : cov::ui::Text::MediumBallStick,
                       language);
}

struct OrbitalAppearanceText {
    const char* material;
    const char* standard;
    const char* glass;
    const char* surface;
    const char* solid;
    const char* wire;
    const char* solid_wire;
    const char* auto_light;
};

OrbitalAppearanceText orbital_appearance_text(const cov::ui::Language language) {
    switch (language) {
        case cov::ui::Language::ChineseSimplified:
            return {"轨道材质", "标准", "玻璃", "表面模式", "实体", "线框", "实体 + 线框", "柔和自动打光"};
        case cov::ui::Language::Japanese:
            return {"軌道マテリアル", "標準", "ガラス", "表示モード", "ソリッド", "ワイヤー", "ソリッド + ワイヤー", "ソフト自動照明"};
        case cov::ui::Language::French:
            return {"Matériau orbital", "Standard", "Verre", "Mode de surface", "Solide", "Filaire", "Solide + filaire", "Éclairage automatique doux"};
        default:
            return {"Orbital material", "Standard", "Glass", "Surface mode", "Solid", "Wire", "Solid + Wire", "Soft automatic lighting"};
    }
}

const char* orbital_material_name(const cov::OrbitalMaterial material,
                                  const OrbitalAppearanceText& text) {
    return material == cov::OrbitalMaterial::Glass ? text.glass : text.standard;
}

const char* orbital_surface_name(const cov::OrbitalSurfaceMode mode,
                                 const OrbitalAppearanceText& text) {
    switch (mode) {
        case cov::OrbitalSurfaceMode::Wire: return text.wire;
        case cov::OrbitalSurfaceMode::SolidWire: return text.solid_wire;
        default: return text.solid;
    }
}

std::optional<std::filesystem::path> file_browser(cov::ui::Language language) {
    std::optional<std::filesystem::path> selected;
    ImGui::SetNextWindowSize(ImVec2(640, 420), ImGuiCond_FirstUseEver);
    if (!ImGui::BeginPopupModal("##file_browser", nullptr, ImGuiWindowFlags_NoSavedSettings)) return selected;
    static auto directory = std::filesystem::current_path();
    static std::array<char, 2048> location{};
    static std::string error;
    ImGui::TextUnformatted(cov::ui::tr(cov::ui::Text::OpenFile, language));
    ImGui::TextWrapped("%s", path_to_utf8(directory).c_str());
    if (ImGui::Button("..") && directory.has_parent_path()) directory = directory.parent_path();
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-1);
    if (ImGui::InputText("##directory", location.data(), location.size(), ImGuiInputTextFlags_EnterReturnsTrue)) {
        std::error_code ec;
        auto next = path_from_utf8(location.data());
        if (std::filesystem::is_directory(next, ec)) { directory = next; error.clear(); }
        else if (std::filesystem::is_regular_file(next, ec)) { selected = next; ImGui::CloseCurrentPopup(); }
        else error = ec ? ec.message() : path_to_utf8(next);
    }
    if (!error.empty()) ImGui::TextWrapped("%s", error.c_str());
    if (ImGui::BeginChild("##entries", ImVec2(0, -36), true)) {
        std::error_code ec;
        std::vector<std::pair<bool, std::filesystem::path>> entries;
        for (std::filesystem::directory_iterator it(directory, std::filesystem::directory_options::skip_permission_denied, ec), end;
             !ec && it != end; it.increment(ec)) {
            const bool folder = it->is_directory(ec);
            if (ec) break;
            entries.emplace_back(folder, it->path());
        }
        if (ec) error = ec.message();
        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
            return a.first != b.first ? a.first > b.first : a.second.filename() < b.second.filename();
        });
        for (const auto& [folder, path] : entries) {
            const auto label = path_to_utf8(path.filename()) + (folder ? "/" : "");
            if (ImGui::Selectable(label.c_str())) {
                if (folder) { directory = path; error.clear(); break; }
                selected = path;
                ImGui::CloseCurrentPopup();
                break;
            }
        }
    }
    ImGui::EndChild();
    const char* close = "Close";
    switch (language) {
        case cov::ui::Language::ChineseSimplified: close = "关闭"; break;
        case cov::ui::Language::Japanese: close = "閉じる"; break;
        case cov::ui::Language::French: close = "Fermer"; break;
        default: break;
    }
    if (ImGui::Button(close)) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
    return selected;
}

} // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
    // Match the UTF-8 path contract used by the file picker and GLFW drops.
    // Narrow CRT argv can otherwise replace non-ANSI characters before parsing.
    std::vector<std::string> utf8_arguments;
    std::vector<char*> argument_pointers;
    int wide_argc = 0;
    if (auto** wide_argv = CommandLineToArgvW(GetCommandLineW(), &wide_argc)) {
        utf8_arguments.reserve(wide_argc);
        for (int i = 0; i < wide_argc; ++i) {
            const int length = WideCharToMultiByte(CP_UTF8, 0, wide_argv[i], -1, nullptr, 0, nullptr, nullptr);
            std::string value(static_cast<std::size_t>(std::max(1, length)), '\0');
            WideCharToMultiByte(CP_UTF8, 0, wide_argv[i], -1, value.data(), length, nullptr, nullptr);
            value.resize(value.size() - 1);
            utf8_arguments.push_back(std::move(value));
        }
        LocalFree(wide_argv);
        for (auto& value : utf8_arguments) argument_pointers.push_back(value.data());
        argument_pointers.push_back(nullptr);
        argc = wide_argc;
        argv = argument_pointers.data();
    }
#endif
    cov::ComputeOptions compute_options;
    std::string input_path;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg.starts_with("--compute-backend=")) compute_options.backend = arg.substr(18);
            else if (arg.starts_with("--compute-device=")) {
                const auto value = arg.substr(17);
                std::size_t end = 0;
                compute_options.device_index = std::stoi(value, &end);
                if (end != value.size() || compute_options.device_index < 0)
                    throw std::invalid_argument("Compute device index must be nonnegative");
            } else if (arg == "--validation-plan" || arg == "--validation-output") ++i;
            else if (!arg.starts_with("--") && input_path.empty()) input_path = arg;
        }
    } catch (const std::exception& e) { std::fprintf(stderr, "%s\n", e.what()); return 2; }
    try { cov::validation::configure(argc, argv); }
    catch (const std::exception& e) { std::fprintf(stderr,"Validation: %s\n",e.what()); return 2; }
    if (!glfwInit()) {
        std::fprintf(stderr, "GLFW initialisation failed\n");
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
    if (cov::validation::background()) {
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        glfwWindowHint(GLFW_FOCUSED, GLFW_FALSE);
    }

    GLFWwindow* window = glfwCreateWindow(
        1500, 940, "Chemical Orbital Visualiser", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "Unable to create OpenGL window\n");
        glfwTerminate();
        return 1;
    }

    glfwSetWindowSizeLimits(window, 640, 360, GLFW_DONT_CARE, GLFW_DONT_CARE);
    glfwMakeContextCurrent(window);
    if (cov::validation::active()) {
        glfwSetWindowSize(window, cov::validation::window_width(), cov::validation::window_height());
        glfwSetWindowTitle(window, "COV native validation");
    }
    glfwSwapInterval(1);
    glfwSetDropCallback(window, drop_callback);

    float x_scale = 1.0f;
    float y_scale = 1.0f;
    glfwGetWindowContentScale(window, &x_scale, &y_scale);
    const float ui_scale = std::clamp(std::max(x_scale, y_scale), 1.0f, 1.75f);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& startup_io = ImGui::GetIO();
    if (cov::validation::active()) {
        startup_io.IniFilename = nullptr;
        // Each native-plan step supplies an ordered event batch for this
        // frame. Do not defer part of it into the next GLFW polling batch.
        startup_io.ConfigInputTrickleEventQueue = false;
    }
    startup_io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    cov::ui::apply_theme(ui_scale);
    cov::ui::configure_fonts(16.5f * ui_scale);
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    int exit_code = 0;
    try {
        cov::VolumeRenderer renderer;
        cov::OrbitCamera camera;

        std::optional<cov::Wavefunction> wavefunction;
        std::optional<cov::OrbitalTrackingResult> frame_tracking;
        std::unique_ptr<cov::OrbitalEvaluator> evaluator;
        cov::GridBox grid_box;

        cov::ui::Language language = cov::ui::Language::English;
        cov::ui::OrbitalUIState orbital_ui;
        cov::MoleculeRenderSettings molecule_render;
        cov::OrbitalMaterial orbital_material = cov::OrbitalMaterial::Standard;
        cov::OrbitalSurfaceMode orbital_surface_mode = cov::OrbitalSurfaceMode::Solid;
        StatusKind status = StatusKind::Ready;
        std::string status_detail;

        std::size_t mo_index = 0;
        std::optional<std::size_t> pending_mo_index;
        float isovalue = 0.03f;
        int resolution = 128;
        std::array<char, 2048> path_buffer{};
        std::filesystem::path current_file;
        std::vector<std::filesystem::path> recent_files;

        bool recompute = false;
        bool resize_and_recompute = false;

        auto evaluate_now = [&]() {
            if (!wavefunction || !evaluator || wavefunction->orbitals.empty()) return;
            if (resize_and_recompute || renderer.nx() != resolution) {
                evaluator->detach_gl_texture();
                renderer.resize_volume(resolution, resolution, resolution);
                evaluator->attach_gl_texture(renderer.volume_texture());
                resize_and_recompute = false;
            }
            if (cov::validation::active()) evaluator->evaluate(mo_index, grid_box, resolution, resolution, resolution);
            else evaluator->begin_evaluate(mo_index, grid_box, resolution, resolution, resolution);
            if (evaluator->ready()) cov::validation::evaluated(mo_index,"selection-or-grid",evaluator->last_kernel_ms());
            status = evaluator->busy() ? StatusKind::Computing : StatusKind::GridUpdated;
            status_detail = evaluator->device_name();
            recompute = false;
        };

        auto load_file = [&](const std::filesystem::path& path) {
            try {
                status = StatusKind::Parsing;
                status_detail = path_to_utf8(path);

                cov::MoldenParseOptions options;
                options.max_atoms = 100;
                options.require_orbitals = true;
                auto wf = cov::parse_molden(path, options);
                if (cov::validation::active()) {
                    std::ostringstream diagnostics;
                    cov::write_numerical_diagnostics_json(diagnostics,wf);
                    cov::validation::record("input.numerical_diagnostics",diagnostics.str());
                    std::ostringstream density_evidence;
                    cov::write_density_evidence_json(density_evidence,wf);
                    cov::validation::record("input.density_evidence",density_evidence.str());
                    std::ostringstream topology_evidence;
                    cov::write_pi_topology_assignments_json(topology_evidence,wf);
                    cov::validation::record("input.pi_topology_evidence",topology_evidence.str());
                }
                const auto new_mo = initial_orbital(wf);
                const auto new_box = make_grid_box(wf);
                std::optional<cov::OrbitalTrackingResult> new_tracking;
                if (wavefunction) {
                    // Cross-frame identity is descriptive state only. Both
                    // canonical wavefunctions remain immutable, and loading a
                    // new frame still resets selection to that frame's own HOMO.
                    new_tracking = cov::track_orbital_subspaces(*wavefunction, wf);
                }

                if (evaluator) evaluator->detach_gl_texture();
                evaluator.reset();
                wavefunction = std::move(wf);
                frame_tracking = std::move(new_tracking);
                renderer.invalidate_geometry_cache();
                evaluator = std::make_unique<cov::OrbitalEvaluator>(*wavefunction, compute_options);
                mo_index = new_mo;
                pending_mo_index.reset();
                orbital_ui.browser_cache={};
                orbital_ui.diagram_cache={};
                grid_box = new_box;

                renderer.resize_volume(resolution, resolution, resolution);
                evaluator->attach_gl_texture(renderer.volume_texture());
                if (cov::validation::active()) evaluator->evaluate(mo_index, grid_box, resolution, resolution, resolution);
                else evaluator->begin_evaluate(mo_index, grid_box, resolution, resolution, resolution);
                if (evaluator->ready()) cov::validation::evaluated(mo_index,"input-load",evaluator->last_kernel_ms());
                current_file = path;
                copy_path_to_buffer(path, path_buffer);
                push_recent(recent_files, path);
                status = evaluator->busy() ? StatusKind::Computing : StatusKind::Loaded;
                status_detail = path_to_utf8(path.filename());
            } catch (const std::exception& e) {
                status = StatusKind::Error;
                status_detail = e.what();
            }
        };

        if (!input_path.empty()) {
            const std::string p = input_path;
            std::snprintf(path_buffer.data(), path_buffer.size(), "%s", p.c_str());
            load_file(path_from_utf8(p));
        }
        if (cov::validation::active() && !wavefunction) {
            throw std::runtime_error("Native validation input failed: "+status_detail);
        }

        bool scene_drag_active = false;

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            if (evaluator) {
                try {
                    if (evaluator->poll()) {
                        status = StatusKind::GridUpdated;
                        status_detail = evaluator->device_name();
                    }
                } catch (const std::exception& e) {
                    status = StatusKind::Error;
                    status_detail = e.what();
                }
            }
            cov::validation::begin_frame(camera,molecule_render,isovalue,resolution,resize_and_recompute);
            if (resize_and_recompute) recompute = true;

            if (!g_dropped_path.empty()) {
                load_file(path_from_utf8(g_dropped_path));
                g_dropped_path.clear();
            }

            if (cov::validation::active()) {
                int window_width=0, window_height=0;
                glfwGetWindowSize(window,&window_width,&window_height);
                if (window_width!=cov::validation::window_width() ||
                    window_height!=cov::validation::window_height()) {
                    glfwSetWindowSize(window,cov::validation::window_width(),cov::validation::window_height());
                    glfwPollEvents();
                }
            }
            int fb_w = 0, fb_h = 0;
            glfwGetFramebufferSize(window, &fb_w, &fb_h);
            if (cov::validation::active() && (fb_w != cov::validation::window_width() ||
                                             fb_h != cov::validation::window_height())) {
                throw std::runtime_error("Validation framebuffer differs from the requested size");
            }
            if (fb_w <= 0 || fb_h <= 0) {
                glfwWaitEventsTimeout(0.05);
                continue;
            }
            ImGui_ImplOpenGL2_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            cov::validation::input_frame();
            ImGui::NewFrame();
            ImGuiIO& io = ImGui::GetIO();
            const auto layout = cov::viewer_layout(io.DisplaySize.x, io.DisplaySize.y,
                                                    fb_w, fb_h, ui_scale);
            const auto& viewport = layout.framebuffer;
            const bool over_scene = layout.scene.contains(io.MousePos.x, io.MousePos.y);
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                scene_drag_active = over_scene && !io.WantCaptureMouse;
            }
            if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) scene_drag_active = false;
            if (scene_drag_active && over_scene && !io.WantCaptureMouse) {
                camera.yaw += io.MouseDelta.x * 0.007f;
                camera.pitch = std::clamp(camera.pitch + io.MouseDelta.y * 0.007f, -1.45f, 1.45f);
            }
            if (over_scene && !io.WantCaptureMouse && std::abs(io.MouseWheel) > 0.0f) {
                camera.distance *= std::pow(0.88f, io.MouseWheel);
                camera.distance = std::clamp(camera.distance, 1.1f, 6.0f);
            }
            glViewport(0, 0, fb_w, fb_h);
            glClearColor(0.025f, 0.031f, 0.043f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glViewport(viewport.x, viewport.y, viewport.width, viewport.height);
            if (wavefunction && viewport.width > 0 && viewport.height > 0) {
                // Geometry and the actual orbital texture share one viewport,
                // projection and depth buffer, outside the control panel.
                renderer.render_geometry(*wavefunction, grid_box, viewport.width, viewport.height,
                                         camera, molecule_render);
                if (evaluator && evaluator->ready()) renderer.render_volume(viewport.width, viewport.height, isovalue, camera,
                                       molecule_render.orbital_opacity,
                                       orbital_material, orbital_surface_mode);
                cov::validation::after_scene(renderer, grid_box, mo_index);
            }
            cov::validation::scene_view(layout, camera);
            glViewport(0, 0, fb_w, fb_h);

            ImGui::SetNextWindowPos(ImVec2(layout.controls.x, layout.controls.y), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(layout.controls.width, layout.controls.height), ImGuiCond_Always);
            ImGui::SetNextWindowBgAlpha(0.965f);
            constexpr ImGuiWindowFlags panel_flags =
                ImGuiWindowFlags_NoTitleBar |
                ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoSavedSettings;

            ImGui::Begin("##cov_control_panel", nullptr, panel_flags);
            const auto panel_position = ImGui::GetWindowPos();
            const auto panel_size = ImGui::GetWindowSize();
            cov::validation::hit("layout.control-panel", panel_position,
                ImVec2(panel_position.x + panel_size.x, panel_position.y + panel_size.y));

            // Keep the complete panel reachable when the window is short.
            ImGui::BeginChild("##cov_panel_scroll", ImVec2(0, 0), false,
                              ImGuiWindowFlags_None);
            const auto brand = [&] {
                ImGui::TextWrapped("%s", cov::ui::tr(cov::ui::Text::AppTitle, language));
                disabled_wrapped(cov::ui::tr(cov::ui::Text::Tagline, language));
            };
            const auto language_control = [&] {
                disabled_wrapped(cov::ui::tr(cov::ui::Text::LanguageLabel, language));
                int language_index = static_cast<int>(language);
                ImGui::SetNextItemWidth(-1.0f);
                if (ImGui::Combo("##language_combo", &language_index,
                                 "English\0简体中文\0日本語\0Français\0")) {
                    language = static_cast<cov::ui::Language>(language_index);
                    glfwSetWindowTitle(window,
                        cov::ui::tr(cov::ui::Text::AppTitle, language));
                }
                cov::validation::item("language");
            };
            const float header_width = ImGui::CalcTextSize(
                cov::ui::tr(cov::ui::Text::AppTitle, language)).x +
                142.0f * ui_scale + 4.0f * ImGui::GetStyle().ItemSpacing.x;
            if (ImGui::GetContentRegionAvail().x < header_width) {
                brand();
                ImGui::Spacing();
                language_control();
            } else if (ImGui::BeginTable("##cov_header", 2,
                                  ImGuiTableFlags_SizingStretchProp |
                                  ImGuiTableFlags_NoSavedSettings)) {
                ImGui::TableSetupColumn("##brand", ImGuiTableColumnFlags_WidthStretch);
                ImGui::TableSetupColumn("##language", ImGuiTableColumnFlags_WidthFixed,
                                        142.0f * ui_scale);
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                brand();
                ImGui::TableNextColumn();
                language_control();
                ImGui::EndTable();
            }

            ImGui::Spacing();
            cov::ui::status_badge(status_label(status, language), status_tone(status));
            if (!status_detail.empty()) {
                disabled_wrapped(status_detail.c_str());
            } else {
                disabled_wrapped(cov::ui::tr(cov::ui::Text::IdleHint, language));
            }
            ImGui::Separator();
            ImGui::Spacing();

            cov::ui::begin_card("##file_card", 192.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::FileSection, language));
            std::optional<std::filesystem::path> recent_to_load;
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::OpenFile, language),
                              ImVec2(150.0f * ui_scale, 0.0f))) {
                const cov::FileDialogResult dialog = cov::open_molden_file_dialog();
                if (dialog.selected()) {
                    load_file(dialog.path);
                } else if (!dialog.supported) {
                    ImGui::OpenPopup("##file_browser");
                } else if (!dialog.cancelled && !dialog.error.empty()) {
                    status = StatusKind::Error;
                    status_detail = dialog.supported
                                        ? dialog.error
                                        : cov::ui::tr(cov::ui::Text::OpenDialogUnsupported, language);
                }
            }
            if (const auto selected = file_browser(language)) load_file(*selected);
            if (!current_file.empty()) {
                const std::string file_label = std::string(
                    cov::ui::tr(cov::ui::Text::CurrentFile, language)) + ": " +
                    path_to_utf8(current_file.filename());
                disabled_wrapped(file_label.c_str());
            }

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoldenPath, language));
            const float load_width = 76.0f * ui_scale;
            ImGui::SetNextItemWidth(std::max(100.0f,
                ImGui::GetContentRegionAvail().x - load_width - 8.0f));
            ImGui::InputText("##molden_path", path_buffer.data(), path_buffer.size());
            ImGui::SameLine();
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::Load, language),
                              ImVec2(load_width, 0.0f))) {
                load_file(path_from_utf8(path_buffer.data()));
            }

            if (!recent_files.empty()) {
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::RecentFiles, language));
                const std::string preview = path_to_utf8(recent_files.front().filename());
                ImGui::SetNextItemWidth(-1.0f);
                if (ImGui::BeginCombo("##recent_files", preview.c_str())) {
                    for (std::size_t i = 0; i < recent_files.size(); ++i) {
                        ImGui::PushID(static_cast<int>(i));
                        const std::string label = path_to_utf8(recent_files[i].filename());
                        if (ImGui::Selectable(label.c_str(), i == 0)) {
                            recent_to_load = recent_files[i];
                        }
                        ImGui::PopID();
                    }
                    ImGui::EndCombo();
                }
            }
            cov::ui::end_card();
            if (recent_to_load) load_file(*recent_to_load);
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##wavefunction_card", 270.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::WavefunctionSection, language));
            if (wavefunction) {
                if (ImGui::BeginTable("##wavefunction_metrics", 2,
                                      ImGuiTableFlags_SizingStretchProp |
                                      ImGuiTableFlags_NoSavedSettings)) {
                    const std::string atoms = std::to_string(wavefunction->atoms.size()) + " / 100";
                    const std::string shells = std::to_string(wavefunction->shells.size());
                    const std::string basis = std::to_string(wavefunction->basis_count);
                    const std::string orbitals = std::to_string(wavefunction->orbitals.size());
                    const std::string convention =
                        std::string("D=") + (wavefunction->pure_d ? "5D" : "6D") +
                        "  F=" + (wavefunction->pure_f ? "7F" : "10F") +
                        "  G=" + (wavefunction->pure_g ? "9G" : "15G");
                    std::string state = "—";
                    if (wavefunction->charge_provenance != cov::DataProvenance::Unavailable ||
                        wavefunction->multiplicity_provenance != cov::DataProvenance::Unavailable) {
                        state.clear();
                        if (wavefunction->charge_provenance != cov::DataProvenance::Unavailable) {
                            if (wavefunction->charge > 0) state += "+";
                            state += std::to_string(wavefunction->charge);
                        } else {
                            state += "?";
                        }
                        state += " / ";
                        state += wavefunction->multiplicity_provenance !=
                                     cov::DataProvenance::Unavailable
                                     ? std::to_string(wavefunction->multiplicity) : "?";
                    }
                    const std::string electron_split =
                        wavefunction->electron_counts_provenance ==
                                cov::DataProvenance::Unavailable
                            ? "— / —"
                            : std::to_string(wavefunction->alpha_electrons) + " / " +
                                  std::to_string(wavefunction->beta_electrons);
                    std::string diagnostics = "—";
                    if (wavefunction->scf_convergence !=
                            cov::ScfConvergenceStatus::Unavailable ||
                        wavefunction->stability !=
                            cov::WavefunctionStabilityStatus::Unavailable) {
                        const char* scf = wavefunction->scf_convergence ==
                                                  cov::ScfConvergenceStatus::Converged
                                              ? cov::ui::tr(cov::ui::Text::Converged,language)
                                              : wavefunction->scf_convergence ==
                                                        cov::ScfConvergenceStatus::Failed
                                                    ? cov::ui::tr(cov::ui::Text::Failed,language)
                                                    : "—";
                        const char* stability = wavefunction->stability ==
                                                        cov::WavefunctionStabilityStatus::Stable
                                                    ? cov::ui::tr(cov::ui::Text::Stable,language)
                                                    : wavefunction->stability ==
                                                              cov::WavefunctionStabilityStatus::Unstable
                                                          ? cov::ui::tr(cov::ui::Text::Unstable,language)
                                                          : "—";
                        diagnostics = std::string(scf) + " / " + stability;
                    }
                    std::string spin_squared = "—";
                    if (wavefunction->spin_squared_provenance !=
                        cov::DataProvenance::Unavailable) {
                        char value[64]{};
                        std::snprintf(value, sizeof(value), "%.4f / %.4f",
                                      wavefunction->spin_squared_before_annihilation,
                                      wavefunction->spin_squared_after_annihilation);
                        spin_squared = value;
                    }
                    metric_row(cov::ui::tr(cov::ui::Text::Atoms, language), atoms.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::Shells, language), shells.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::BasisFunctions, language), basis.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::Orbitals, language), orbitals.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::ShellConvention, language), convention.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::ChargeMultiplicity, language),
                               state.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::AlphaBetaElectrons, language),
                               electron_split.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::SCFStability, language),
                               diagnostics.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::SpinSquared, language),
                               spin_squared.c_str());
                    ImGui::EndTable();
                }
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##frame_tracking_card", 156.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::FrameTracking,
                                               language));
            if (frame_tracking) {
                if (ImGui::BeginTable("##frame_tracking_metrics", 2,
                                      ImGuiTableFlags_SizingStretchProp |
                                      ImGuiTableFlags_NoSavedSettings)) {
                    std::size_t matched_members = 0u;
                    for (const auto& match : frame_tracking->matches) {
                        matched_members += match.from_members.size();
                    }
                    const std::string matched =
                        std::to_string(frame_tracking->matches.size()) +
                        " (" + std::to_string(matched_members) + " MO)";
                    const std::string unmatched =
                        std::to_string(frame_tracking->unmatched_from.size()) +
                        " / " +
                        std::to_string(frame_tracking->unmatched_to.size());
                    metric_row(cov::ui::tr(
                                   cov::ui::Text::AtomMappingCompatibility,
                                   language),
                               cov::ui::tr(
                                   frame_tracking->atom_mapping_compatible
                                       ? cov::ui::Text::Compatible
                                       : cov::ui::Text::Incompatible,
                                   language));
                    metric_row(cov::ui::tr(cov::ui::Text::MatchedSubspaces,
                                           language),
                               matched.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::UnmatchedSubspaces,
                                           language),
                               unmatched.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::TrackingOptimisation,
                                           language),
                               cov::ui::tr(
                                   frame_tracking->composite_optimisation_truncated
                                       ? cov::ui::Text::ConservativeFallback
                                       : cov::ui::Text::ExactOrNotNeeded,
                                   language));
                    ImGui::EndTable();
                }
            } else {
                ImGui::TextDisabled("%s", cov::ui::tr(
                    cov::ui::Text::NoPreviousFrame, language));
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::validation::anchor("panel.browser");
            cov::ui::begin_card("##orbital_browser_card", 620.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::OrbitalBrowser, language));
            cov::ui::OrbitalUIActions orbital_actions;
            if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                cov::ui::draw_orbital_browser(*wavefunction, mo_index, orbital_ui,
                                              language, ui_scale, orbital_actions);
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::validation::anchor("panel.diagram");
            cov::ui::begin_card("##energy_diagram_card", 430.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::EnergyDiagram, language));
            cov::ui::OrbitalUIActions diagram_actions;
            if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                cov::ui::draw_energy_diagram(*wavefunction, mo_index, orbital_ui,
                                             language, ui_scale, diagram_actions);
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            if (orbital_actions.select_orbital) pending_mo_index = orbital_actions.select_orbital;
            if (diagram_actions.select_orbital) pending_mo_index = diagram_actions.select_orbital;
            const bool export_requested = orbital_actions.export_diagram || diagram_actions.export_diagram;
            if (export_requested && wavefunction) {
                std::filesystem::path base = current_file.empty()
                                                 ? std::filesystem::current_path() / "mo_diagram"
                                                 : current_file;
                base = cov::validation::export_base(base);
                const auto snapshot=diagram_actions.drawn_diagram;
                cov::MODiagramExportResult result;
                if (snapshot) result=cov::export_mo_diagram_bundle(*snapshot,base);
                else result.error="No current diagram view is available for export";
#ifdef COV_ENABLE_VALIDATION
                cov::validation::record("export.actual","{\"base\":"+cov::validation::quote(path_to_utf8(base))+
                    ",\"snapshot_id\":"+(snapshot?cov::validation::quote(snapshot->data.view->id):"null")+
                    ",\"mode\":"+(snapshot?std::to_string(static_cast<int>(snapshot->data.mode)):"null")+
                    ",\"selected_index\":"+(snapshot && snapshot->data.view->inspected_orbital_index
                        ?std::to_string(*snapshot->data.view->inspected_orbital_index):"null")+
                    ",\"success\":"+((result.svg&&result.png&&result.json&&result.csv)?"true":"false")+"}");
#endif
                if (result.svg && result.png && result.json && result.csv) {
                    status = StatusKind::Exported;
                    status_detail = path_to_utf8(result.svg_path.parent_path() /
                        result.svg_path.stem()) + ".{png,svg,json,csv}";
                } else {
                    status = StatusKind::Error;
                    status_detail = result.error.empty()
                                        ? cov::ui::tr(cov::ui::Text::ExportFailed, language)
                                        : result.error;
                }
            }

            cov::ui::begin_card("##render_card", 545.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::RenderingSection, language));
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoleculeStyle, language));
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::BeginCombo("##molecule_style",
                                  molecule_style_name(molecule_render.style, language))) {
                for (const cov::MoleculeStyle style : {
                         cov::MoleculeStyle::MediumBallAndStick,
                         cov::MoleculeStyle::StickDelocalisation}) {
                    const bool selected = molecule_render.style == style;
                    if (ImGui::Selectable(molecule_style_name(style, language), selected)) {
                        molecule_render.style = style;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            if (molecule_render.style == cov::MoleculeStyle::StickDelocalisation) {
                ImGui::TextDisabled("%s",
                    cov::ui::tr(cov::ui::Text::DelocalisationHeuristic, language));
            }

            const OrbitalAppearanceText appearance = orbital_appearance_text(language);
            ImGui::TextDisabled("%s", appearance.material);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::BeginCombo("##orbital_material",
                                  orbital_material_name(orbital_material, appearance))) {
                for (const cov::OrbitalMaterial material : {
                         cov::OrbitalMaterial::Standard,
                         cov::OrbitalMaterial::Glass}) {
                    const bool selected = orbital_material == material;
                    if (ImGui::Selectable(orbital_material_name(material, appearance), selected)) {
                        orbital_material = material;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }

            ImGui::TextDisabled("%s", appearance.surface);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::BeginCombo("##orbital_surface",
                                  orbital_surface_name(orbital_surface_mode, appearance))) {
                for (const cov::OrbitalSurfaceMode mode : {
                         cov::OrbitalSurfaceMode::Solid,
                         cov::OrbitalSurfaceMode::Wire,
                         cov::OrbitalSurfaceMode::SolidWire}) {
                    const bool selected = orbital_surface_mode == mode;
                    if (ImGui::Selectable(orbital_surface_name(mode, appearance), selected)) {
                        orbital_surface_mode = mode;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            ImGui::TextDisabled("%s", appearance.auto_light);

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::AtomSize, language));
            ImGui::SliderFloat("##atom_size", &molecule_render.atom_scale, 0.55f, 1.8f, "%.2f");
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::BondSize, language));
            ImGui::SliderFloat("##bond_size", &molecule_render.bond_scale, 0.5f, 2.0f, "%.2f");
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowHydrogens, language),
                            &molecule_render.show_hydrogens);
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowCoordinationContacts, language),
                            &molecule_render.show_coordination_contacts);
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowMulticentreSupport, language),
                            &molecule_render.show_multicentre_support);
            ImGui::Checkbox(cov::ui::tr(
                                cov::ui::Text::ShowPolyhedralCageSupport,language),
                            &molecule_render.show_polyhedral_cage_support);
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowWeakInteractions, language),
                            &molecule_render.show_weak_interactions);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s", cov::ui::tr(cov::ui::Text::WeakInteractionsHint,
                                                     language));
            }

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoleculeOpacity, language));
            ImGui::SliderFloat("##molecule_opacity", &molecule_render.molecule_opacity,
                               0.15f, 1.0f, "%.2f");
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::OrbitalOpacity, language));
            ImGui::SliderFloat("##orbital_opacity", &molecule_render.orbital_opacity,
                               0.02f, 1.0f, "%.2f");

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::Isovalue, language));
            ImGui::SliderFloat("##isovalue", &isovalue, 0.002f, 0.12f, "%.4f",
                               ImGuiSliderFlags_Logarithmic);

            constexpr int resolutions[] = {64, 128, 256, 512};
            int resolution_index = 1;
            for (int i = 0; i < 4; ++i) {
                if (resolutions[i] == resolution) resolution_index = i;
            }
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::Grid, language));
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::Combo("##grid", &resolution_index,
                             "64³\0" "128³\0" "256³\0" "512³\0")) {
                resolution = resolutions[resolution_index];
                resize_and_recompute = true;
                recompute = true;
            }

            const float button_gap = ImGui::GetStyle().ItemSpacing.x;
            const float half_button = (ImGui::GetContentRegionAvail().x - button_gap) * 0.5f;
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::RecomputeGrid, language),
                              ImVec2(half_button, 0.0f))) {
                recompute = true;
            }
            ImGui::SameLine();
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::ResetCamera, language),
                              ImVec2(half_button, 0.0f))) {
                camera = {};
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##performance_card", 190.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::PerformanceSection, language));
            if (evaluator) {
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::CUDADevice, language));
                ImGui::TextUnformatted(evaluator->device_name());
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::LastKernel, language));
                ImGui::Text("%.3f ms", evaluator->last_kernel_ms());
                if (evaluator->ready()) cov::ui::status_badge(cov::ui::tr(cov::ui::Text::GPUResident, language), cov::ui::Tone::Success);
                else if (evaluator->busy()) disabled_wrapped(status_label(StatusKind::Computing, language));
            } else {
                ImGui::TextDisabled("%s —", cov::ui::tr(cov::ui::Text::CUDADevice, language));
            }
            ImGui::TextDisabled("%s: %s",
                                cov::ui::tr(cov::ui::Text::FontStatus, language),
                                cov::ui::font_status());
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::InteractionHint, language));
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::IsovalueHint, language));
            cov::ui::end_card();

            ImGui::EndChild();
            ImGui::End();
            cov::validation::field("language",std::to_string(static_cast<int>(language)));
            cov::validation::ui_frame(mo_index,pending_mo_index.value_or(mo_index));

            // Selection debounce: at most the latest requested orbital is evaluated
            // once at the end of this frame. Browser hover/filtering does not recompute the grid.
            if (pending_mo_index && wavefunction &&
                *pending_mo_index < wavefunction->orbitals.size()) {
                if (*pending_mo_index != mo_index) {
                    mo_index = *pending_mo_index;
                    recompute = true;
                }
                pending_mo_index.reset();
            }

            if (recompute) {
                try {
                    evaluate_now();
                } catch (const std::exception& e) {
                    status = StatusKind::Error;
                    status_detail = e.what();
                    recompute = false;
                }
            }

            ImGui::Render();
            ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
            cov::validation::end_frame(fb_w,fb_h,mo_index,orbital_ui,wavefunction?&*wavefunction:nullptr);

            glfwSwapBuffers(window);
            if (cov::validation::done()) {exit_code=cov::validation::result();break;}
        }

        if (evaluator) evaluator->detach_gl_texture();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Fatal error: %s\n", e.what());
        exit_code = 1;
    }

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return exit_code;
}
