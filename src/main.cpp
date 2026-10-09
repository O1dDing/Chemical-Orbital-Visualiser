#include "cov/orbital_evaluator.hpp"
#include "cov/numerical_diagnostics.hpp"
#include "cov/pi_topology_evidence.hpp"
#include <sstream>
#include "cov/file_dialog.hpp"
#include "cov/gl_api.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/molden_parser.hpp"
#include "cov/molecule_style.hpp"
#include "cov/nbo_ui.hpp"
#include "cov/nbo_molecular_overlay.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/orbital_tracking.hpp"
#include "cov/nbo_aomo_labels.hpp"
#include "cov/orbital_symmetry_components.hpp"
#include "cov/nbo_aomo_text.hpp"
#include "cov/orbital_ui.hpp"
#include "cov/ui.hpp"
#include "cov/ui_raster_text.hpp"
#include "cov/volume_renderer.hpp"
#include "cov/validation.hpp"
#include "cov/viewer_layout.hpp"
#include "cov/ui_window_layers.hpp"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <shellapi.h>
#endif

namespace {

std::vector<std::filesystem::path> g_dropped_paths;

void report_fatal_error(const char* message) {
    std::fprintf(stderr,"COV: %s\n",message);
#ifdef _WIN32
    if(!cov::validation::active()) {
        const int count=MultiByteToWideChar(CP_UTF8,0,message,-1,nullptr,0);
        std::wstring text(static_cast<std::size_t>(std::max(1,count)),L'\0');
        MultiByteToWideChar(CP_UTF8,0,message,-1,text.data(),count);
        MessageBoxW(nullptr,text.c_str(),L"Chemical Orbital Visualiser",MB_OK|MB_ICONERROR);
    }
#endif
}

enum class StatusKind {
    Ready,
    Parsing,
    Loaded,
    GridUpdated,
    Exported,
    Error,
    Computing,
};

void drop_callback(GLFWwindow*, const int count, const char** paths) {
    if (!paths) return;
    g_dropped_paths.clear();
    for (int i=0;i<count;++i) if(paths[i])
        g_dropped_paths.push_back(std::filesystem::u8path(paths[i]));
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

const char* scene_text(cov::ui::Language language, const char* en, const char* zh,
                      const char* ja, const char* fr) {
    switch(language) {
        case cov::ui::Language::ChineseSimplified: return zh;
        case cov::ui::Language::Japanese: return ja;
        case cov::ui::Language::French: return fr;
        default: return en;
    }
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
        report_fatal_error("GLFW initialisation failed");
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
        report_fatal_error("Unable to create OpenGL window");
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

        std::unique_ptr<cov::Wavefunction> wavefunction;
        std::optional<cov::Wavefunction> nbo_wavefunction;
        std::optional<cov::NboIntegration> integration;
        std::optional<cov::RoutedAnalysis> routed;
        std::unique_ptr<cov::NboSelectionView> inspection;
        cov::InteractionGraph semantic_graph;
        // Each overlay retains a separate signed field and compute evaluator.
        struct AdditionalField {
            std::unique_ptr<cov::VolumeRenderer> renderer;
            std::unique_ptr<cov::OrbitalEvaluator> evaluator;
            std::size_t index=0;
        };
        std::vector<AdditionalField> additional_fields;
        cov::ui::NboUIState nbo_ui;
        bool nbo_active = false;
        std::size_t canonical_mo_index = 0;
        std::size_t nbo_mo_index = 0;
        std::optional<cov::OrbitalTrackingResult> frame_tracking;
        std::unique_ptr<cov::OrbitalEvaluator> evaluator;
        cov::GridBox grid_box;

        cov::ui::Language language = cov::ui::Language::English;
        cov::ui::OrbitalUIState orbital_ui;
        orbital_ui.nbo_ui = &nbo_ui;
        cov::MoleculeRenderSettings molecule_render;
        cov::OrbitalMaterial orbital_material = cov::OrbitalMaterial::Standard;
        cov::OrbitalSurfaceMode orbital_surface_mode = cov::OrbitalSurfaceMode::Solid;
        StatusKind status = StatusKind::Ready;
        std::string status_detail;
        std::string status_error_detail;
        const bool profile_open=std::getenv("COV_PROFILE_OPEN")!=nullptr;
        bool profile_frame_pending=false;
        auto profile_start=std::chrono::steady_clock::now();
        auto profile_previous=profile_start;
        const auto profile_stage=[&](const char* stage) {
            if(!profile_open)return;
            const auto now=std::chrono::steady_clock::now();
            std::fprintf(stderr,"COV open: %s stage=%.3f total=%.3f seconds\n",stage,
                std::chrono::duration<double>(now-profile_previous).count(),
                std::chrono::duration<double>(now-profile_start).count());
            std::fflush(stderr);profile_previous=now;
        };

        std::size_t mo_index = 0;
        std::optional<std::size_t> pending_mo_index;
        float isovalue = 0.03f;
        int resolution = 128;
        std::array<char, 2048> path_buffer{};
        std::filesystem::path current_file;
        std::vector<std::filesystem::path> recent_files;

        bool recompute = false;
        bool resize_and_recompute = false;

        auto start_grid = [&](cov::OrbitalEvaluator& target, std::size_t index,
                              const cov::GridBox& box) {
            if (cov::validation::active())
                target.evaluate(index, box, resolution, resolution, resolution);
            else
                target.begin_evaluate(index, box, resolution, resolution, resolution);
        };
        auto fields_busy = [&]() {
            return (evaluator && evaluator->busy()) ||
                std::any_of(additional_fields.begin(), additional_fields.end(),
                    [](const auto& field) { return field.evaluator->busy(); });
        };
        auto fields_ready = [&]() {
            return evaluator && evaluator->ready() &&
                std::all_of(additional_fields.begin(), additional_fields.end(),
                    [](const auto& field) { return field.evaluator->ready(); });
        };

        auto active_wavefunction = [&]() -> const cov::Wavefunction* {
            if(inspection) return &inspection->wavefunction;
            return nbo_active && nbo_wavefunction ? &*nbo_wavefunction
                                                   : wavefunction ? &*wavefunction : nullptr;
        };
        auto active_view = [&]() {
            cov::ActiveOrbitalView view;
            if(inspection) {
                view.kind=cov::ActiveOrbitalKind::Inspection;
                view.source_id=inspection->selection.dataset_id;
                view.label=inspection->label;
                view.source_label=inspection->selection.label;
                view.display_name_evidence="verified selected combination and source terms";
                view.semantic_kind=inspection->selection.semantic_kind;
                view.group_id=inspection->selection.group_id;
                view.selection=inspection->selection;
                view.canonical_index=inspection->selection.target_canonical_index;
                view.rendered_index=mo_index;
                if(!inspection->selection.terms.empty())
                    view.spin=inspection->selection.terms.front().orbital.spin;
                view.source_spin=view.spin;
                view.spin_semantics="actual selected source orbital spin";
                if(inspection->selection.spatial_spin){
                    view.spin=cov::NboSpin::Total;
                    view.spin_semantics="verified common spatial orbital; both source spins retained in selection; source_spin identifies the rendering basis";
                    view.display_name_evidence=inspection->selection.spatial_spin->correspondence;
                    if(inspection->selection.semantic_kind=="salc_component"&&view.canonical_index&&wavefunction&&*view.canonical_index<wavefunction->orbitals.size()){
                        view.spin=wavefunction->orbitals[*view.canonical_index].spin==cov::Spin::Beta?cov::NboSpin::Beta:cov::NboSpin::Alpha;
                        view.spin_semantics="target canonical component spin; source_spin identifies the verified common spatial rendering basis";
                    }
                }
                const auto& selected=inspection->selection;
                if(selected.mode==cov::NboSelectionMode::Orbital && selected.terms.size()==1 &&
                   selected.terms.front().orbital.kind==cov::NboOrbitalKind::Canonical && wavefunction) {
                    const auto index=selected.terms.front().orbital.index;
                    const auto* name=nbo_ui.aomo.names && index<nbo_ui.aomo.names->canonical.size()?
                        &nbo_ui.aomo.names->canonical[index]:nullptr;
                    view.label=cov::ui::canonical_mo_display_label(*wavefunction,index,name);
                    if(name){view.display_name_metadata_json=cov::ui::serialize_orbital_name_json(*name);if(name->verified)view.display_name_evidence=name->detail;}
                } else if(selected.semantic_kind=="canonical_symmetry_component" &&
                          selected.target_canonical_index && wavefunction) {
                    const auto names=nbo_ui.aomo.names?nbo_ui.aomo.names:cov::ui::canonical_mo_names(*wavefunction);
                    const auto index=*selected.target_canonical_index;
                    if(names && index<names->canonical.size()) {
                        view.display_name_metadata_json=cov::ui::serialize_orbital_name_json(names->canonical[index]);
                        view.display_name_evidence=names->canonical[index].decomposition_status;
                    }
                } else if(nbo_ui.aomo.salc_model && nbo_ui.aomo.names &&
                          selected.dataset_id==nbo_ui.aomo.salc_model->dataset_id &&
                          !selected.source_id.empty()) {
                    const auto& model=*nbo_ui.aomo.salc_model;
                    for(std::size_t i=0;i<model.orbitals.size() && i<nbo_ui.aomo.names->salc.size();++i)
                        if(model.orbitals[i].id==selected.source_id) {
                            const auto& name=nbo_ui.aomo.names->salc[i];
                            view.display_name_metadata_json=cov::ui::serialize_orbital_name_json(name);
                            if(selected.semantic_kind=="salc_symmetry_component") {
                                view.spin=model.orbitals[i].spin;
                                view.spin_semantics="source SALC spin; source_spin and signed terms identify the verified spatial reconstruction basis";
                                view.display_name_evidence=name.decomposition_status;
                            }
                            if(selected.semantic_kind=="salc" || selected.semantic_kind=="spin_averaged_spatial_orbital") {
                                if(nbo_ui.aomo.drawn_snapshot && nbo_ui.aomo.drawn_snapshot->integration_id==selected.dataset_id)
                                    for(const auto& node:nbo_ui.aomo.drawn_snapshot->nodes)
                                        if(node.id=="salc:"+selected.source_id) {
                                            view.label=node.individual_label.empty()?node.label:node.individual_label;
                                            if(!node.name_detail.empty())view.display_name_evidence=node.name_detail;
                                            break;
                                        }
                            }
                            if(name.verified && !name.irrep.empty()) {
                                // The actual drawn object is the title. Atomic side orbitals
                                // keep their atom/shell name; symmetry is a separate field.
                                if(selected.semantic_kind=="salc_component")
                                    view.label=name.label+scene_text(language," · component"," · 分量"," · 成分"," · composante");
                                view.display_name_evidence=name.detail;
                            }
                            break;
                        }
                }
            }else if(nbo_active && nbo_ui.dataset &&
                     mo_index<nbo_ui.dataset->orbitals.size()) {
                view.kind=cov::ActiveOrbitalKind::NboSet;
                const auto& row=nbo_ui.dataset->orbitals[mo_index];
                view.source_id=integration?integration->id:nbo_ui.dataset->source.path;
                view.label=row.label;
                view.source_label=row.label;
                view.display_name_evidence="associated producer NBO orbital";
                view.semantic_kind="nbo";
                view.spin=row.spin;
                view.source_spin=row.spin;
                view.spin_semantics="producer NBO spin block";
                view.group_id="nbo:"+std::string(cov::nbo_spin_name(row.spin))+":"+
                    std::to_string(row.id);
                if(integration && row.id) {
                    const cov::NboOrbitalRef ref{cov::NboOrbitalKind::NBO,row.spin,row.id-1};
                    if(cov::nbo_orbital(*integration,ref))
                        view.selection=cov::nbo_single_selection(*integration,ref);
                }
                view.rendered_index=mo_index;
            }else if(wavefunction && mo_index<wavefunction->orbitals.size()) {
                view.kind=cov::ActiveOrbitalKind::Canonical;
                view.source_id=routed?routed->canonical_fingerprint:path_to_utf8(current_file);
                view.source_label=cov::ui::canonical_mo_source_label(*wavefunction,mo_index);
                const cov::ui::NboAomoName* display_name=nullptr;
                const auto current_names=nbo_ui.aomo.names?nbo_ui.aomo.names:cov::ui::canonical_mo_names(*wavefunction);
                if(current_names && mo_index<current_names->canonical.size()){
                    const auto& name=current_names->canonical[mo_index];
                    display_name=&name;
                    view.display_name_metadata_json=cov::ui::serialize_orbital_name_json(name);
                    if(name.verified)view.display_name_evidence=name.detail;
                }
                view.label=cov::ui::canonical_mo_display_label(*wavefunction,mo_index,display_name);
                view.semantic_kind="canonical";
                view.spin=wavefunction->orbital_occupation_model==cov::OrbitalOccupationModel::CanonicalShared?
                    cov::NboSpin::Total:wavefunction->orbitals[mo_index].spin==cov::Spin::Beta?
                    cov::NboSpin::Beta:cov::NboSpin::Alpha;
                view.source_spin=wavefunction->orbitals[mo_index].spin==cov::Spin::Beta?
                    cov::NboSpin::Beta:cov::NboSpin::Alpha;
                view.spin_semantics=view.spin==cov::NboSpin::Total?
                    "shared spatial-orbital display; alpha is the producer channel field":
                    "explicit canonical spin orbital";
                view.canonical_index=mo_index;
                view.rendered_index=mo_index;
            }
            return view;
        };
        auto export_analysis_companions = [&](std::filesystem::path base) {
            if(!routed || !wavefunction)return;
            base.replace_extension();
            const auto write=[&](const char* suffix,const std::string& payload) {
                auto path=base;path+=suffix;
                std::ofstream out(path,std::ios::binary);
                if(!out)throw std::runtime_error("Cannot write routed analysis companion");
                out<<payload;
                if(!out)throw std::runtime_error("Routed analysis companion write failed");
            };
            std::vector<std::size_t> export_members;
            if(nbo_ui.aomo.drawn_snapshot)export_members=nbo_ui.aomo.drawn_snapshot->central_mo_indices;
            else if(orbital_ui.diagram_cache.snapshot) {
                for(const auto& level:orbital_ui.diagram_cache.snapshot->data.levels) {
                    export_members.insert(export_members.end(),level.member_indices.begin(),level.member_indices.end());
                    export_members.insert(export_members.end(),level.member_spin_counterparts.begin(),level.member_spin_counterparts.end());
                }
            } else if(const auto selected=active_view().canonical_index)export_members.push_back(*selected);
            const auto analysis=cov::serialize_routed_analysis_json(*routed,false,false,export_members);
            const auto active=cov::serialize_active_orbital_view_json(active_view());
            write(".analysis.json",analysis);
            write(".active-view.json",active);
            // The current diagram JSON already owns its names, source
            // mappings and display values. Do not attach whole input matrices
            // and full-model duplicate files to a current-view export.
            cov::validation::record("chemistry.export",
                "{\"analysis\":"+analysis+",\"active_view\":"+active+"}");
        };

        auto identity = [&]() {
            const auto* active = active_wavefunction();
            const auto* orbital = active && mo_index < active->orbitals.size()
                                    ? &active->orbitals[mo_index] : nullptr;
            const auto spin = inspection && !inspection->selection.terms.empty() ? inspection->selection.terms.front().orbital.spin :
                              nbo_active && nbo_ui.dataset &&
                              mo_index<nbo_ui.dataset->orbitals.size()
                ? nbo_ui.dataset->orbitals[mo_index].spin
                : orbital && orbital->spin == cov::Spin::Beta
                    ? cov::NboSpin::Beta : cov::NboSpin::Alpha;
            const cov::NboCanonicalEvidence* evidence=nullptr;
            if(nbo_ui.dataset)for(const auto& record:nbo_ui.dataset->association.canonical_evidence)
                if(record.spin==spin || (spin==cov::NboSpin::Alpha && record.spin==cov::NboSpin::Total)){
                    evidence=&record;break;
                }
            const bool single_inspection=inspection && inspection->selection.terms.size()==1;
            const auto* descriptor=single_inspection && integration?
                cov::nbo_orbital(*integration,inspection->selection.terms[0].orbital):nullptr;
            const bool canonical_inspection=inspection && inspection->selection.mode==cov::NboSelectionMode::Orbital &&
                descriptor && descriptor->ref.kind==cov::NboOrbitalKind::Canonical;
            const bool direct_canonical=wavefunction && wavefunction->source==cov::WavefunctionSource::Fchk &&
                (!inspection || canonical_inspection) && !nbo_active;
            cov::validation::orbital_identity(
                inspection ? (inspection->selection.mode==cov::NboSelectionMode::Orbital && descriptor?
                    cov::nbo_orbital_kind_name(descriptor->ref.kind):"signed-combination") : nbo_active ? "nbo" : "canonical",
                inspection && integration ? integration->id : nbo_active && nbo_ui.dataset ? nbo_ui.dataset->source.path
                                             : path_to_utf8(current_file),
                cov::nbo_spin_name(spin),
                descriptor ? descriptor->ref.index : inspection ? std::numeric_limits<std::size_t>::max() : orbital ? orbital->source_orbital_index
                        : std::numeric_limits<std::size_t>::max(),
                nbo_ui.dataset ? nbo_ui.dataset->association.status : "not_attached",
                descriptor ? descriptor->source.path : inspection ? "verified signed orbital terms" :
                    direct_canonical ? path_to_utf8(current_file) : evidence ? evidence->coefficient_source : "not_available",
                direct_canonical,
                evidence && evidence->density_verified);
        };

        auto clear_inspection_controls = [&] {
            nbo_ui.inspected_nlmo.reset();nbo_ui.inspected_nho.reset();
            nbo_ui.nho_sum_owner.reset();nbo_ui.nho_sum_nao_indices.clear();
            nbo_ui.inspected_dataset_id.clear();
        };
        auto activate_set = [&](bool use_nbo) {
            if (use_nbo == nbo_active && !inspection) {clear_inspection_controls();return;}
            if (use_nbo && !nbo_wavefunction) throw std::runtime_error("NBO coefficients are not available for rendering");
            const cov::Wavefunction& target = use_nbo ? *nbo_wavefunction : *wavefunction;
            if (target.orbitals.empty()) throw std::runtime_error("Selected orbital set is empty");
            auto next_evaluator = std::make_unique<cov::OrbitalEvaluator>(target, compute_options);
            const std::size_t next_index = std::min(use_nbo ? nbo_mo_index : canonical_mo_index,
                                                    target.orbitals.size()-1);
            if (evaluator) evaluator->detach_gl_texture();
            const bool resized=resize_and_recompute || renderer.nx()!=resolution;
            try {
                if(resized)renderer.resize_volume(resolution,resolution,resolution);
                next_evaluator->attach_gl_texture(renderer.volume_texture());
                start_grid(*next_evaluator, next_index, grid_box);
            } catch (...) {
                next_evaluator->detach_gl_texture();
                if (evaluator) {
                    evaluator->attach_gl_texture(renderer.volume_texture());
                    start_grid(*evaluator, mo_index, grid_box);
                }
                throw;
            }
            evaluator = std::move(next_evaluator);
            additional_fields.clear();
            inspection.reset();
            nbo_ui.aomo.selection.reset();
            nbo_ui.aomo.pending_selection.reset();
            nbo_ui.aomo.selected_side_node_id.clear();
            ++nbo_ui.aomo.revision;
            clear_inspection_controls();
            nbo_active = use_nbo;
            mo_index = next_index;
            pending_mo_index.reset();
            orbital_ui.browser_cache={};
            orbital_ui.diagram_cache={};
            renderer.invalidate_geometry_cache();
            resize_and_recompute=false;
            identity();
            if (fields_ready()) cov::validation::evaluated(mo_index,"set-switch",evaluator->last_kernel_ms());
            status=fields_busy() ? StatusKind::Computing : StatusKind::GridUpdated;
            status_error_detail.clear();
            status_detail=use_nbo ? "NBO "+std::to_string(next_index+1) :
                cov::ui::canonical_mo_display_label(*wavefunction,next_index,
                    nbo_ui.aomo.names && next_index<nbo_ui.aomo.names->canonical.size()
                        ? &nbo_ui.aomo.names->canonical[next_index] : nullptr);
            recompute = false;
        };

        auto evaluate_now = [&]() {
            const auto* active = active_wavefunction();
            if (!active || !evaluator || active->orbitals.empty()) return;
            if (resize_and_recompute || renderer.nx() != resolution) {
                evaluator->detach_gl_texture();
                renderer.resize_volume(resolution, resolution, resolution);
                evaluator->attach_gl_texture(renderer.volume_texture());
                resize_and_recompute = false;
            }
            start_grid(*evaluator, mo_index, grid_box);
            for(auto& field:additional_fields){
                if(field.renderer->nx()!=resolution){
                    field.evaluator->detach_gl_texture();
                    field.renderer->resize_volume(resolution,resolution,resolution);
                    field.evaluator->attach_gl_texture(field.renderer->volume_texture());
                }
                start_grid(*field.evaluator, field.index, grid_box);
            }
            identity();
            if (fields_ready()) cov::validation::evaluated(mo_index,"selection-or-grid",evaluator->last_kernel_ms());
            status = fields_busy() ? StatusKind::Computing : StatusKind::GridUpdated;
            status_error_detail.clear();
            status_detail = evaluator->device_name();
            recompute = false;
        };

        auto load_file = [&](const std::filesystem::path& path) {
            try {
                status = StatusKind::Parsing;
                status_detail = path_to_utf8(path);
                status_error_detail.clear();

                cov::MoldenParseOptions options;
                options.max_atoms = 100;
                options.require_orbitals = true;
                auto next_wavefunction=std::make_unique<cov::Wavefunction>(cov::parse_molden(path, options));
                profile_stage("wavefunction");
                auto& wf=*next_wavefunction;
                if (cov::validation::active() && !cov::validation::forensic_mode()) {
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
                    // TODO(perf) ✳: Profile tracking separately from parsing and
                    // first display when investigating slow subsequent file opens.
                    // Cross-frame identity is descriptive state only. Both
                    // canonical wavefunctions remain immutable, and loading a
                    // new frame still resets selection to that frame's own HOMO.
                    // Track with a shared work/deadline budget; never retain a
                    // partial assignment as a completed correspondence.
                    new_tracking = cov::track_orbital_subspaces(*wavefunction, wf);
                    profile_stage("frame-matching");
                    if(cov::validation::active()){
                        std::ostringstream tracking;
                        tracking<<"{\"budget_exhausted\":"<<(new_tracking->tracking_budget_exhausted?"true":"false")
                            <<",\"budget_reason\":"<<int(new_tracking->budget_exhaustion_reason)
                            <<",\"budget_stage\":"<<int(new_tracking->budget_exhausted_stage)
                            <<",\"work_units\":"<<new_tracking->tracking_work_units
                            <<",\"matches\":"<<new_tracking->matches.size()
                            <<",\"unmatched_from\":"<<new_tracking->unmatched_from.size()
                            <<",\"unmatched_to\":"<<new_tracking->unmatched_to.size()
                            <<",\"unresolved_from\":"<<new_tracking->unresolved_from.size()
                            <<",\"unresolved_to\":"<<new_tracking->unresolved_to.size()<<'}';
                        cov::validation::record("input.frame_tracking",tracking.str());
                    }
                }

                auto next_route=cov::route_chemistry(wf);
                auto next_graph=*next_route.interaction_graph.value;
                profile_stage("orbital-analysis");
                auto next_evaluator=std::make_unique<cov::OrbitalEvaluator>(wf, compute_options);
                if (evaluator) evaluator->detach_gl_texture();
                try {
                    renderer.resize_volume(resolution,resolution,resolution);
                    next_evaluator->attach_gl_texture(renderer.volume_texture());
                    start_grid(*next_evaluator, new_mo, new_box);
                } catch(...) {
                    next_evaluator->detach_gl_texture();
                    if(evaluator){evaluator->attach_gl_texture(renderer.volume_texture());
                        start_grid(*evaluator, mo_index, grid_box);}
                    throw;
                }
                evaluator=std::move(next_evaluator);
                profile_stage("orbital-grid");
                additional_fields.clear();
                inspection.reset();
                nbo_ui.integration=nullptr;
                integration.reset();
                routed.reset();
                nbo_ui.aomo={};clear_inspection_controls();
                nbo_ui.selected_atoms.clear();
                nbo_ui.selected_structure.reset();
                nbo_ui.atom_colour_mode=0;
                nbo_ui.show_bond_indices=false;nbo_ui.show_e2=false;
                nbo_ui.input_discovery.reset();nbo_ui.pending_candidate.reset();
                nbo_active = false;
                nbo_wavefunction.reset();
                nbo_ui.dataset.reset();
                nbo_ui.focus = {};
                nbo_ui.error.clear();
                nbo_ui.export_status.clear();
                wavefunction = std::move(next_wavefunction);
                cov::ui::invalidate_canonical_mo_names_cache();
                routed=std::move(next_route);
                nbo_ui.routed=&*routed;
                semantic_graph=std::move(next_graph);
                frame_tracking = std::move(new_tracking);
                renderer.invalidate_geometry_cache();
                mo_index = new_mo;
                canonical_mo_index = new_mo;
                nbo_mo_index = 0;
                pending_mo_index.reset();
                orbital_ui.browser_cache={};
                orbital_ui.diagram_cache={};
                if(cov::validation::active() && !cov::validation::forensic_mode())
                    cov::validation::record("chemistry.route",cov::serialize_routed_analysis_json(*routed));
                grid_box = new_box;
                current_file = path;

                identity();
                if (fields_ready()) cov::validation::evaluated(mo_index,"input-load",evaluator->last_kernel_ms());
                copy_path_to_buffer(path, path_buffer);
                push_recent(recent_files, path);
                status = evaluator->busy() ? StatusKind::Computing : StatusKind::Loaded;
                status_detail = path_to_utf8(path.filename());
            } catch (const std::exception& e) {
                status = StatusKind::Error;
                status_error_detail=e.what();
                status_detail=scene_text(language,"The calculation file could not be read.",
                    "无法读取计算文件。","計算ファイルを読み込めません。",
                    "Impossible de lire le fichier de calcul.");
            }
        };

        auto apply_selection = [&](const cov::NboOrbitalSelection& selection) {
            if(!wavefunction) throw std::runtime_error("No canonical wavefunction is loaded");
            std::optional<cov::NboIntegration> canonical_components;
            const cov::NboIntegration* source=integration?&*integration:nullptr;
            if(!source && selection.semantic_kind=="canonical_symmetry_component") {
                canonical_components=cov::ui::canonical_component_dataset(*wavefunction);
                source=&*canonical_components;
            }
            if(!source)throw std::runtime_error("No verified orbital data is attached");
            auto next=std::make_unique<cov::NboSelectionView>(cov::make_nbo_selection_view(*source,*wavefunction,selection));
            if(!next->available || next->wavefunction.orbitals.empty())
                throw std::runtime_error(next->detail.empty()?next->status:next->detail);
            auto next_evaluator=std::make_unique<cov::OrbitalEvaluator>(next->wavefunction, compute_options);
            std::vector<AdditionalField> extra;
            for(std::size_t i=1;i<next->wavefunction.orbitals.size();++i){
                AdditionalField f;
                f.renderer=std::make_unique<cov::VolumeRenderer>();
                f.renderer->resize_volume(resolution,resolution,resolution);
                f.evaluator=std::make_unique<cov::OrbitalEvaluator>(next->wavefunction, compute_options);
                f.evaluator->attach_gl_texture(f.renderer->volume_texture());
                start_grid(*f.evaluator, i, grid_box);
                f.index=i;extra.push_back(std::move(f));
            }
            if(evaluator)evaluator->detach_gl_texture();
            try {
                if(renderer.nx()!=resolution)renderer.resize_volume(resolution,resolution,resolution);
                next_evaluator->attach_gl_texture(renderer.volume_texture());
                start_grid(*next_evaluator, 0, grid_box);
            } catch(...) {
                next_evaluator->detach_gl_texture();
                if(evaluator){evaluator->attach_gl_texture(renderer.volume_texture());
                    start_grid(*evaluator, mo_index, grid_box);}
                throw;
            }
            evaluator=std::move(next_evaluator);
            additional_fields=std::move(extra);
            inspection=std::move(next);
            nbo_active=false;mo_index=0;
            nbo_ui.aomo.selection=selection;
            // Keep parent controls only for that actual parent/component inspection.
            const auto is_kind=[&](cov::NboOrbitalKind kind){return selection.mode==cov::NboSelectionMode::Orbital &&
                selection.terms.size()==1 && selection.terms.front().orbital.kind==kind;};
            if(is_kind(cov::NboOrbitalKind::NLMO))nbo_ui.inspected_nlmo=selection.terms.front().orbital;
            else if(selection.group_id.rfind("nlmo:",0)!=0)nbo_ui.inspected_nlmo.reset();
            if(is_kind(cov::NboOrbitalKind::NHO))nbo_ui.inspected_nho=selection.terms.front().orbital;
            else if(selection.group_id.rfind("nho:",0)!=0){
                nbo_ui.inspected_nho.reset();nbo_ui.nho_sum_owner.reset();nbo_ui.nho_sum_nao_indices.clear();
            }
            nbo_ui.aomo.selected_side_node_id.clear();
            if(nbo_ui.aomo.drawn_snapshot)for(const auto& node:nbo_ui.aomo.drawn_snapshot->nodes){
                const bool salc=!selection.source_id.empty() && node.salc_index &&
                    node.id=="salc:"+selection.source_id;
                const bool local=selection.terms.size()==1 && node.orbital &&
                    node.lane!=cov::ui::NboAomoLane::Centre && *node.orbital==selection.terms.front().orbital;
                if(salc || local){nbo_ui.aomo.selected_side_node_id=node.id;break;}
            }
            ++nbo_ui.aomo.revision;
            status_detail=inspection->label;
            // Match the scene's canonical identity rather than an ambiguous
            // spin-local descriptor number (for example beta MO 3 vs MO 31).
            if(selection.terms.size()==1 &&
               selection.terms.front().orbital.kind==cov::NboOrbitalKind::Canonical) {
                const auto& ref=selection.terms.front().orbital;
                status_detail=cov::ui::canonical_mo_display_label(*wavefunction,ref.index,
                    nbo_ui.aomo.names && ref.index<nbo_ui.aomo.names->canonical.size()
                        ? &nbo_ui.aomo.names->canonical[ref.index] : nullptr);
            }
            nbo_ui.aomo.status=status_detail;
            nbo_ui.error.clear();
            status=fields_busy() ? StatusKind::Computing : StatusKind::GridUpdated;
            status_error_detail.clear();
            nbo_ui.selected_atoms={inspection->atoms.begin(),inspection->atoms.end()};
            if(selection.target_canonical_index)canonical_mo_index=*selection.target_canonical_index;
            else if(selection.terms.size()==1 && selection.terms.front().orbital.kind==cov::NboOrbitalKind::Canonical)
                canonical_mo_index=selection.terms.front().orbital.index;
            pending_mo_index.reset();recompute=false;resize_and_recompute=false;
            identity();
            if(cov::validation::active() && !cov::validation::forensic_mode())
                cov::validation::record("aomo.selection",cov::serialize_nbo_selection_json(*inspection));
            if (fields_ready()) cov::validation::evaluated(0,"typed-orbital-selection",evaluator->last_kernel_ms());
        };

        auto attach_integration = [&](cov::NboIntegration next,bool keep_discovery=false) {
            if(!wavefunction)throw std::runtime_error("Load the matching Gaussian wavefunction first");
            cov::annotate_nbo_bond_channels(next,*wavefunction);
            if(inspection || nbo_active)activate_set(false);
            nbo_wavefunction.reset();
            integration=std::move(next);
            nbo_ui.integration=&*integration;
            nbo_ui.dataset=integration->dataset;
            routed=cov::route_chemistry(*wavefunction,&*integration);
            nbo_ui.routed=&*routed;
            semantic_graph=*routed->interaction_graph.value;
            renderer.invalidate_geometry_cache();
            if(cov::validation::active() && !cov::validation::forensic_mode())
                cov::validation::record("chemistry.route",cov::serialize_routed_analysis_json(*routed));
            nbo_ui.focus={};nbo_ui.aomo={};clear_inspection_controls();
            nbo_ui.selected_atoms.clear();nbo_ui.selected_structure.reset();
            nbo_ui.show_bond_indices=false;nbo_ui.show_e2=false;
            nbo_ui.pending_candidate.reset();
            if(!keep_discovery)nbo_ui.input_discovery.reset();
            nbo_ui.input_status=integration->dataset.source.path;
            nbo_ui.error.clear();
            if(const auto* c=cov::nbo_capability(*integration,"nbo");c && c->available()) {
                try { nbo_wavefunction=cov::make_nbo_wavefunction(integration->dataset,*wavefunction); }
                catch(const std::exception& e){nbo_ui.error=e.what();}
            }
            const auto suggested=cov::ui::suggest_initial_nbo_index(*integration,orbital_ui.diagram_cache.snapshot.get());
            nbo_ui.selected_orbital=suggested.value_or(std::numeric_limits<std::size_t>::max());
            nbo_mo_index=nbo_ui.selected_orbital;
            // Prepare the attachment's single name set before the scene asks
            // for its caption. Otherwise that first caption builds a second,
            // standalone canonical name set immediately before the attached one.
            cov::ui::prepare_nbo_aomo_state(nbo_ui.aomo,*integration,*wavefunction);
            status=nbo_ui.error.empty()?StatusKind::Loaded:StatusKind::Error;
            status_detail=nbo_ui.error.empty()?integration->dataset.source.path:nbo_ui.error;
            // Do not build large diagnostic payloads during ordinary viewing.
            if(cov::validation::active() && !cov::validation::forensic_mode())
                cov::validation::record("nbo.integration",cov::serialize_nbo_integration_json(*integration));
            cov::validation::record("nbo.attach","{\"source\":"+cov::validation::quote(integration->dataset.source.path)+
                ",\"association\":"+cov::validation::quote(integration->dataset.association.status)+
                ",\"renderable\":"+(nbo_wavefunction?"true":"false")+"}");
        };
        auto clear_integration = [&] {
            if(inspection || nbo_active)activate_set(false);
            nbo_wavefunction.reset();nbo_ui.integration=nullptr;integration.reset();
            nbo_ui.dataset.reset();nbo_ui.aomo={};nbo_ui.focus={};clear_inspection_controls();
            nbo_ui.selected_atoms.clear();nbo_ui.selected_structure.reset();
            nbo_ui.selected_orbital=std::numeric_limits<std::size_t>::max();
            nbo_ui.atom_colour_mode=0;nbo_ui.show_bond_indices=false;nbo_ui.show_e2=false;
            nbo_ui.pending_candidate.reset();
            orbital_ui.browser_cache={};orbital_ui.diagram_cache={};
            routed.reset();nbo_ui.routed=nullptr;semantic_graph={};
            if(wavefunction){
                routed=cov::route_chemistry(*wavefunction);
                nbo_ui.routed=&*routed;semantic_graph=*routed->interaction_graph.value;
                if(cov::validation::active() && !cov::validation::forensic_mode())
                    cov::validation::record("chemistry.route",cov::serialize_routed_analysis_json(*routed));
            }
            renderer.invalidate_geometry_cache();
        };
        auto apply_candidate = [&](const cov::NboInputCandidate& candidate) {
            if(!candidate.canonical.empty() && candidate.canonical.lexically_normal()!=current_file.lexically_normal()) {
                load_file(candidate.canonical);
                if(status==StatusKind::Error)throw std::runtime_error(status_detail);
            }
            if(!wavefunction)throw std::runtime_error("The package has no loaded canonical wavefunction");
            if(candidate.report.empty()){
                clear_integration();
                nbo_ui.input_status=scene_text(language,"NBO data could not be located.",
                    "无法定位 NBO 数据。","NBO データが見つかりません。","Données NBO introuvables.");
                status=StatusKind::Loaded;
                status_detail=nbo_ui.input_status;
                return;
            }
            auto next_integration=cov::read_nbo_integration(*wavefunction,candidate);
            profile_stage("nbo-read-and-association");
            attach_integration(std::move(next_integration),true);
            profile_stage("nbo-analysis");
            nbo_ui.input_status=candidate.label;
            if(const auto* matched=cov::nbo_capability(*integration,"source_association");
               matched && matched->state==cov::NboCapabilityState::Rejected){
                nbo_ui.error="NBO data mismatch: "+matched->detail;
                nbo_ui.input_status=scene_text(language,"The NBO data does not match this calculation.",
                    "NBO 数据与当前计算不匹配。","NBO データが現在の計算と一致しません。",
                    "Les données NBO ne correspondent pas à ce calcul.");
                status_error_detail=matched->detail;
                status=StatusKind::Error;status_detail=nbo_ui.input_status;
            }
        };
        auto load_inputs = [&](const std::vector<std::filesystem::path>& paths) {
            profile_start=profile_previous=std::chrono::steady_clock::now();
            profile_frame_pending=true;
            try {
                nbo_ui.error.clear();
                status_error_detail.clear();
                auto found=cov::discover_nbo_inputs(paths);
                profile_stage("input-discovery");
                nbo_ui.pending_candidate.reset();
                if(found.candidates.size()==1 && !found.selection_required) {
                    apply_candidate(found.candidates.front());
                } else if(found.candidates.empty()) {
                    // Existing Molden and standalone Gaussian use the same loader.
                    if(paths.size()==1 && !std::filesystem::is_directory(paths.front())){
                        load_file(paths.front());
                        if(status==StatusKind::Error)throw std::runtime_error(status_detail);
                    }
                    nbo_ui.input_status=scene_text(language,"NBO data could not be located.",
                        "无法定位 NBO 数据。","NBO データが見つかりません。","Données NBO introuvables.");
                } else {
                    nbo_ui.input_status=scene_text(language,"Multiple calculations found.",
                        "找到多个计算。","複数の計算が見つかりました。","Plusieurs calculs trouvés.");
                }
                nbo_ui.input_discovery=std::move(found);
            } catch(const std::exception& e){
                clear_integration();
                nbo_ui.input_discovery.reset();
                nbo_ui.error=e.what();nbo_ui.input_status=scene_text(language,"NBO data could not be loaded.",
                    "无法载入 NBO 数据。","NBO データを読み込めません。","Impossible de charger les données NBO.");
                status=StatusKind::Error;status_detail=nbo_ui.input_status;
                cov::validation::record("input.package.error","{\"reason\":"+cov::validation::quote(e.what())+"}");
            }
        };

        if (!input_path.empty()) {
            const std::string p = input_path;
            std::snprintf(path_buffer.data(), path_buffer.size(), "%s", p.c_str());
            load_inputs({path_from_utf8(p)});
        }
        if (cov::validation::active() && !wavefunction) {
            throw std::runtime_error("Native validation input failed: "+status_detail);
        }

        bool scene_drag_active = false;
        ImVec2 scene_press{};
        bool scene_was_dragged=false;
        bool choose_nbo_input=false;
        unsigned profile_draw_frames=0;

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            // Open at the next frame boundary, after previous draw references
            // are released. Reuse the package loader and native picker.
            if(choose_nbo_input){
                choose_nbo_input=false;
                const auto dialog=cov::open_wavefunction_file_dialog(language,true);
                if(dialog.selected())load_inputs({dialog.path});
                else if(!dialog.cancelled&&!dialog.error.empty()){
                    status=StatusKind::Error;
                    status_detail=dialog.supported?dialog.error:
                        cov::ui::tr(cov::ui::Text::OpenDialogUnsupported,language);
                }
            }
            if (evaluator) {
                try {
                    bool completed = evaluator->poll();
                    for (auto& field : additional_fields)
                        completed = field.evaluator->poll() || completed;
                    if (completed && fields_ready() && status == StatusKind::Computing)
                        status = StatusKind::GridUpdated;
                } catch (const std::exception& e) {
                    // A selected overlay is one view: do not leave a partial
                    // set of fields visible when one compute request fails.
                    evaluator->cancel();
                    for (auto& field : additional_fields) field.evaluator->cancel();
                    status = StatusKind::Error;
                    status_error_detail = e.what();
                    status_detail = scene_text(language,"The orbital grid could not be calculated.",
                        "无法计算轨道网格。","軌道グリッドを計算できませんでした。",
                        "Impossible de calculer la grille orbitale.");
                }
            }
            cov::validation::begin_frame(camera,molecule_render,isovalue,resolution,resize_and_recompute);
            if (resize_and_recompute) recompute = true;
            if(auto paths=cov::validation::take_dropped_paths();!paths.empty())load_inputs(paths);

            if (!g_dropped_paths.empty()) {
                const auto paths=std::move(g_dropped_paths);
                g_dropped_paths.clear();
                load_inputs(paths);
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
            bool canonical_requested=false;
            const auto* aomo_cap=integration?cov::nbo_capability(*integration,"aomo"):nullptr;
            const bool aomo_available=aomo_cap && aomo_cap->available();
            const auto layout = cov::viewer_layout(io.DisplaySize.x, io.DisplaySize.y,
                                                    fb_w, fb_h, ui_scale,aomo_available);
            const auto& viewport = layout.framebuffer;
            const bool over_scene = layout.scene.contains(io.MousePos.x, io.MousePos.y);
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                scene_drag_active = over_scene && !io.WantCaptureMouse;
                scene_press=io.MousePos;scene_was_dragged=false;
            }
            if(scene_drag_active && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                if(!scene_was_dragged && over_scene && integration){
                    const auto picked=renderer.pick_geometry((io.MousePos.x-layout.scene.x)/layout.scene.width,
                        (io.MousePos.y-layout.scene.y)/layout.scene.height);
                    if(picked){
                        if(picked->kind==cov::GeometryTargetKind::Atom){
                            nbo_ui.selected_structure.reset();
                            if(!io.KeyCtrl)nbo_ui.selected_atoms.clear();
                            if(nbo_ui.selected_atoms.count(picked->index))nbo_ui.selected_atoms.erase(picked->index);
                            else nbo_ui.selected_atoms.insert(picked->index);
                        }else if(picked->index<integration->structure.size()){
                            nbo_ui.selected_structure=picked->index;
                            const auto& e=integration->structure[picked->index];
                            nbo_ui.selected_atoms={e.atoms.begin(),e.atoms.end()};
                            if(e.kind=="donor_acceptor" && e.orbitals.size()>=2){
                                cov::NboOrbitalSelection s;s.dataset_id=integration->id;s.label=e.label;
                                s.mode=cov::NboSelectionMode::Overlay;
                                s.terms={{e.orbitals[0],1},{e.orbitals[1],1}};nbo_ui.aomo.pending_selection=s;
                            }
                        }
                        cov::validation::record("scene.pick","{\"kind\":"+
                            std::to_string(static_cast<int>(picked->kind))+",\"index\":"+std::to_string(picked->index)+"}");
                    }
                }
                scene_drag_active=false;
            }
            if(scene_drag_active && std::hypot(io.MousePos.x-scene_press.x,io.MousePos.y-scene_press.y)>4*ui_scale)
                scene_was_dragged=true;
            if (scene_drag_active && scene_was_dragged && ImGui::IsMouseDown(ImGuiMouseButton_Left) && over_scene && !io.WantCaptureMouse) {
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
            identity();
            if(cov::validation::active())
                cov::validation::field("scene.active_view",
                    cov::serialize_active_orbital_view_json(active_view()));
            if(routed) {
                cov::validation::field("chemistry.route.charge",
                    std::string(cov::routed_status_name(routed->total_atomic_charge.status))+":"+
                    cov::routed_provider_name(routed->total_atomic_charge.provider)+":"+
                    routed->total_atomic_charge.fallback_reason);
                cov::validation::field("chemistry.route.spin",
                    std::string(cov::routed_status_name(routed->atomic_spin.status))+":"+
                    cov::routed_provider_name(routed->atomic_spin.provider)+":"+
                    routed->atomic_spin.reason);
            }
            std::optional<cov::MoleculeOverlay> overlay;
            if(integration && wavefunction)overlay=cov::make_nbo_molecule_overlay(*integration,semantic_graph,
                wavefunction->atoms.size(),{nbo_ui.selected_atoms.begin(),nbo_ui.selected_atoms.end()},
                nbo_ui.selected_structure,static_cast<cov::AtomScalarMode>(nbo_ui.atom_colour_mode),
                nbo_ui.show_bond_indices,nbo_ui.show_e2,routed?&*routed:nullptr,
                wavefunction.get(),nbo_ui.show_lewis_skeleton?cov::NboBondDisplayMode::LewisStructure:
                    cov::NboBondDisplayMode::DefaultSkeleton,
                nbo_ui.aomo.salc_model?&nbo_ui.aomo.salc_model->symmetry_scope:nullptr);
            if(overlay && cov::validation::active())cov::validation::field("overlay.scalars",
                cov::serialize_molecule_overlay_scalars_json(*overlay,routed?&*routed:nullptr));
            if (const auto* active = active_wavefunction();
                active && viewport.width > 0 && viewport.height > 0) {
                // Geometry and the actual orbital texture share one viewport,
                // projection and depth buffer, outside the control panel.
                renderer.render_geometry(*wavefunction, grid_box, viewport.width, viewport.height,
                                         camera, molecule_render,overlay?&*overlay:nullptr,
                                         &semantic_graph);
                if (fields_ready()) {
                    renderer.render_volume(viewport.width, viewport.height, isovalue, camera,
                                           molecule_render.orbital_opacity,
                                           orbital_material, orbital_surface_mode);
                    for(auto& field:additional_fields)field.renderer->render_volume(viewport.width,viewport.height,
                        isovalue,camera,molecule_render.orbital_opacity,orbital_material,orbital_surface_mode,1);
                    cov::validation::after_scene(renderer, grid_box, mo_index);
                    for(std::size_t f=0;f<additional_fields.size();++f)
                        cov::validation::after_scene(*additional_fields[f].renderer,grid_box,mo_index,f+1);
                }
            }
            cov::validation::scene_view(layout, camera);
            glViewport(0, 0, fb_w, fb_h);

            // A non-intercepting scene legend also supplies validation hit
            // locations. Both ordinary clicks and test clicks use pick_geometry.
            ImGui::SetNextWindowPos(ImVec2(layout.scene.x,layout.scene.y));
            ImGui::SetNextWindowSize(ImVec2(layout.scene.width,layout.scene.height));
            ImGui::SetNextWindowBgAlpha(0);
            ImGui::Begin("##scene_evidence",nullptr,ImGuiWindowFlags_NoInputs|ImGuiWindowFlags_NoDecoration|
                ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoBringToFrontOnFocus);
            if(wavefunction){
                const auto spin_caption=[&](cov::NboSpin spin)->const char* {
                    return spin==cov::NboSpin::Alpha?"α":spin==cov::NboSpin::Beta?"β":
                        scene_text(language,"total","总","全体","total");
                };
                const auto mo_caption=[&](std::size_t index) {
                    const auto* name=nbo_ui.aomo.names && index<nbo_ui.aomo.names->canonical.size()
                        ?&nbo_ui.aomo.names->canonical[index]:nullptr;
                    return cov::ui::canonical_mo_display_label(*wavefunction,index,name);
                };
                std::string label;
                if(inspection){
                    label=inspection->label;
                    if(inspection->selection.terms.size()==1){const auto& ref=inspection->selection.terms[0].orbital;
                        label+=" | "+std::string(cov::nbo_orbital_kind_name(ref.kind))+" "+std::to_string(ref.index+1)+" / "+spin_caption(ref.spin);
                    }else {
                        label+=" | "+std::to_string(inspection->selection.terms.size())+
                            scene_text(language," signed terms"," 个带系数项"," 個の符号付き成分"," termes signés");
                        std::vector<std::string> source_kinds;
                        for(const auto& term:inspection->selection.terms){
                            const auto identity=std::string(cov::nbo_orbital_kind_name(term.orbital.kind))+
                                " / "+spin_caption(term.orbital.spin);
                            if(std::find(source_kinds.begin(),source_kinds.end(),identity)==source_kinds.end())
                                source_kinds.push_back(identity);
                        }
                        for(const auto& identity:source_kinds)label+=" | "+identity;
                    }
                }
                else if(nbo_active && nbo_ui.dataset && mo_index<nbo_ui.dataset->orbitals.size())
                    label="NBO "+std::to_string(mo_index+1)+" / "+spin_caption(nbo_ui.dataset->orbitals[mo_index].spin);
                else label=mo_caption(canonical_mo_index);
                // Translate presentation only; inspection identities and the actual signed field stay intact.
                if(inspection) {
                    const auto& selected=inspection->selection;
                    if(selected.spatial_spin || cov::ui::is_symmetry_component(selected))
                        label=selected.label;
                    else if(selected.mode==cov::NboSelectionMode::Orbital && selected.terms.size()==1 &&
                       selected.terms[0].orbital.kind==cov::NboOrbitalKind::Canonical)
                        label=mo_caption(selected.terms[0].orbital.index);
                    else if(selected.target_canonical_index)
                        label=std::string(scene_text(language,"Orbital components of ","轨道组成：","軌道成分：","Composantes orbitales de "))+
                            mo_caption(*selected.target_canonical_index);
                    else if(selected.mode==cov::NboSelectionMode::Combination && nbo_ui.aomo.salc_model &&
                            selected.dataset_id==nbo_ui.aomo.salc_model->dataset_id) {
                        const auto& model=*nbo_ui.aomo.salc_model;
                        for(std::size_t i=0;i<model.orbitals.size();++i) {
                            const auto& side=model.orbitals[i];
                            if(side.terms.size()!=selected.terms.size())continue;
                            bool same=true;
                            for(std::size_t j=0;j<side.terms.size();++j)
                                same=same && side.terms[j].orbital==selected.terms[j].orbital &&
                                    side.terms[j].coefficient==selected.terms[j].coefficient;
                            if(!same)continue;
                            label=side.atoms.size()==1?"NAO":
                                ((nbo_ui.aomo.names && i<nbo_ui.aomo.names->salc.size())?
                                    nbo_ui.aomo.names->salc[i].label+" · SALC":"SALC");
                            label+=" (";
                            for(std::size_t j=0;j<side.atoms.size();++j) {
                                if(j)label+=", ";const auto atom=side.atoms[j];
                                if(atom<wavefunction->atoms.size())label+=wavefunction->atoms[atom].symbol+std::to_string(atom+1);
                            }
                            label+=") / ";label+=spin_caption(side.spin);break;
                        }
                    }
                }
                if(inspection && (inspection->selection.semantic_kind=="salc" ||
                    inspection->selection.semantic_kind=="spin_averaged_spatial_orbital"))label=active_view().label;
                ImGui::TextWrapped("%s",label.c_str());
                if(inspection && inspection->selection.mode==cov::NboSelectionMode::Overlay){
                    for(std::size_t i=0;i<inspection->selection.terms.size();++i){
                        const auto& term=inspection->selection.terms[i];
                        ImGui::TextWrapped("%zu: %s %zu / %s; c=%+.6g",i+1,
                            cov::nbo_orbital_kind_name(term.orbital.kind),term.orbital.index+1,
                            cov::nbo_spin_name(term.orbital.spin),term.coefficient);
                    }
                }
                ImGui::TextColored(ImVec4(.96f,.5f,.45f,1),"%s",scene_text(language,
                    "Phase: red + / blue -","相位：红 + / 蓝 −","位相：赤 + / 青 −","Phase : rouge + / bleu -"));
                if(!additional_fields.empty())ImGui::TextColored(ImVec4(.3f,.92f,.75f,1),"%s",
                    scene_text(language,"Overlay: green + / gold - (separate field)",
                        "叠加轨道：绿 + / 金 −（独立场）","重ね合わせ：緑 + / 金 −（独立場）",
                        "Superposition : vert + / or - (champ distinct)"));
                if(overlay && std::any_of(overlay->bonds.begin(),overlay->bonds.end(),[](const auto& bond){
                    return bond.style==cov::OverlayBondStyle::Unresolved;}))
                    ImGui::TextWrapped("%s",scene_text(language,
                        "Grey dotted links: connectivity without an assigned integer bond order.",
                        "灰色点划线：已识别的连接，整数键级未确定。",
                        "灰色の点線：整数結合次数が未確定の結合。",
                        "Pointillés gris : connexion sans ordre de liaison entier attribué."));
                if(overlay && overlay->colour_mode!=cov::AtomScalarMode::Element){
                    ImGui::Text("%s: -%.3f ... 0 ... +%.3f",
                        overlay->colour_mode==cov::AtomScalarMode::NaturalCharge?
                            scene_text(language,"NPA charge / e","NPA 电荷 / e","NPA 電荷 / e","Charge NPA / e"):
                            scene_text(language,"Spin population / e","自旋布居 / e","スピン分布 / e","Population de spin / e"),overlay->scalar_range,overlay->scalar_range);
                    ImGui::TextDisabled("%s",scene_text(language,
                        "Atoms: purple (-), grey (0), gold (+); missing stays uncoloured",
                        "原子色标：紫（负）→灰（零）→金（正）；缺失不着色",
                        "原子の色：紫（負）→灰（零）→金（正）；欠損値は着色しません",
                        "Atomes : violet (-), gris (0), or (+) ; sans couleur si absent"));
                }
                for(const auto& t:renderer.geometry_targets()){
                    const float px=layout.scene.x+(t.segment?(t.x+t.end_x)/2:t.x)*layout.scene.width;
                    const float py=layout.scene.y+(t.segment?(t.y+t.end_y)/2:t.y)*layout.scene.height;
                    const std::string kind=t.kind==cov::GeometryTargetKind::Atom?"atom":
                        t.kind==cov::GeometryTargetKind::Bond?"bond":t.kind==cov::GeometryTargetKind::Relation?"relation":"multicentre";
                    cov::validation::hit("scene."+kind+"."+std::to_string(t.index),ImVec2(px-3,py-3),ImVec2(px+3,py+3));
                    if(molecule_render.show_numbers && molecule_render.number_atoms &&
                       t.kind==cov::GeometryTargetKind::Atom && t.index<wavefunction->atoms.size() &&
                       !(molecule_render.number_ignore_h && wavefunction->atoms[t.index].atomic_number==1)) {
                        const auto atom_label=wavefunction->atoms[t.index].symbol+std::to_string(t.index+1);
                        ImGui::GetWindowDrawList()->AddText(ImVec2(px+8,py-18),IM_COL32(224,231,239,255),atom_label.c_str());
                        cov::validation::field("scene.number.atom."+std::to_string(t.index),atom_label);
                    }
                    if(overlay && nbo_ui.show_bond_indices && t.kind==cov::GeometryTargetKind::Bond &&
                       integration && t.index<integration->structure.size()){
                        const auto& e=integration->structure[t.index];
                        if(e.wiberg){char value[64];std::snprintf(value,sizeof(value),"WBI %.3f",*e.wiberg);
                            ImGui::GetWindowDrawList()->AddText(ImVec2(px+8,py+8),IM_COL32(255,225,145,255),value);}
                    }
                }
                if(molecule_render.show_numbers && molecule_render.number_fragments)
                    for(const auto& group:nbo_ui.aomo.fragment_groups) {
                        float x=0,y=0;std::size_t count=0;
                        for(const auto& target:renderer.geometry_targets())
                            if(target.kind==cov::GeometryTargetKind::Atom && group.atoms.contains(target.index) &&
                               target.index<wavefunction->atoms.size() &&
                               !(molecule_render.number_ignore_h && wavefunction->atoms[target.index].atomic_number==1)) {
                                x+=target.x;y+=target.y;++count;
                            }
                        if(!count)continue;
                        const auto text="L"+std::to_string(group.id);
                        ImGui::GetWindowDrawList()->AddText(ImVec2(layout.scene.x+x/count*layout.scene.width,
                            layout.scene.y+y/count*layout.scene.height+16),IM_COL32(255,222,147,255),text.c_str());
                        cov::validation::field("scene.number.fragment."+std::to_string(group.id),text);
                    }
                if(integration && nbo_ui.selected_structure && *nbo_ui.selected_structure<integration->structure.size()){
                    const auto& evidence=integration->structure[*nbo_ui.selected_structure];
                    ImGui::TextWrapped("%s",cov::ui::nbo_structure_display_label(
                        evidence,wavefunction.get(),language).c_str());
                    if(evidence.value)ImGui::Text("%.5g %s",*evidence.value,evidence.units.c_str());
                }
                if(overlay)for(auto atom:nbo_ui.selected_atoms)if(atom<overlay->atom_values.size() && overlay->atom_values[atom])
                    ImGui::Text("Atom %zu: %+.6f",atom+1,*overlay->atom_values[atom]);
                cov::validation::field("scene.orbital_label",label);
                cov::validation::field("scene.overlay_fields",std::to_string(additional_fields.size()+1));
            }
            ImGui::End();

            if(wavefunction && layout.scene.width>220 && layout.scene.height>160){
                const float toolbar_width=std::min(layout.scene.width-24.0f,520.0f*ui_scale);
                const char* opacity_help=scene_text(language,"Orbital opacity (lower to see atoms and bonds)",
                    "轨道不透明度（降低可看清原子与键）","軌道の不透明度（下げると原子と結合が見えます）",
                    "Opacité orbitale (réduire pour voir atomes et liaisons)");
                const char* threshold_help=scene_text(language,
                    "Display threshold only; orbital coefficients and amplitudes stay unchanged.",
                    "仅调整显示阈值；轨道系数与幅度不变。",
                    "表示しきい値のみ調整します。軌道係数と振幅は変わりません。",
                    "Seuil d’affichage seul ; coefficients et amplitudes restent inchangés.");
                const auto& toolbar_style=ImGui::GetStyle();
                const float text_width=std::max(1.0f,toolbar_width-2*toolbar_style.WindowPadding.x);
                const char* reveal_label=scene_text(language,"Reveal bonds","看清骨架","骨格を表示","Voir les liaisons");
                const char* iso_label=scene_text(language,"Isovalue","等值面","等値面","Isovaleur");
                const char* fit_label=scene_text(language,"Fit component","适合当前成分","成分に合わせる","Adapter");
                const float fit_width=ImGui::CalcTextSize(fit_label).x+2*toolbar_style.FramePadding.x;
                const float reveal_width=ImGui::CalcTextSize(reveal_label).x+2*toolbar_style.FramePadding.x;
                const float value_width=ImGui::CalcTextSize("1.2345e-07").x+2*toolbar_style.FramePadding.x;
                const bool stacked=text_width<ImGui::CalcTextSize(iso_label).x+value_width+
                    fit_width+2*toolbar_style.ItemSpacing.x;
                const float toolbar_height=ImGui::CalcTextSize(opacity_help,nullptr,false,text_width).y+
                    ImGui::CalcTextSize(threshold_help,nullptr,false,text_width).y+
                    (stacked?4:2)*ImGui::GetFrameHeight()+(stacked?6:4)*toolbar_style.ItemSpacing.y+
                    (molecule_render.show_numbers?4:1)*ImGui::GetFrameHeightWithSpacing()+
                    2*toolbar_style.WindowPadding.y;
                ImGui::SetNextWindowPos(ImVec2(layout.scene.x+12.0f,
                    std::max(layout.scene.y+12.0f,layout.scene.y+layout.scene.height-toolbar_height-12.0f)));
                ImGui::SetNextWindowSize(ImVec2(toolbar_width,toolbar_height));
                ImGui::SetNextWindowBgAlpha(.82f);
                ImGui::Begin("##scene_display_controls",nullptr,ImGuiWindowFlags_NoDecoration|
                    ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoSavedSettings|
                    cov::ui::background_panel_flags);
                ImGui::TextWrapped("%s",opacity_help);
                ImGui::SetNextItemWidth(stacked?text_width:std::max(70.0f,text_width-reveal_width-toolbar_style.ItemSpacing.x));
                ImGui::SliderFloat("##scene_orbital_opacity",&molecule_render.orbital_opacity,.02f,1.0f,"%.2f");
                cov::validation::item("scene.opacity");
                if(!stacked)ImGui::SameLine();
                if(ImGui::Button(reveal_label))molecule_render.orbital_opacity=.24f;
                cov::validation::item("scene.reveal_bonds");
                ImGui::TextUnformatted(iso_label);ImGui::SameLine();
                ImGui::SetNextItemWidth(std::max(65.0f,ImGui::GetContentRegionAvail().x-
                    (stacked?0:fit_width+toolbar_style.ItemSpacing.x)));
                ImGui::SliderFloat("##scene_isovalue",&isovalue,1e-7f,.2f,"%.5g",ImGuiSliderFlags_Logarithmic);
                cov::validation::item("scene.isovalue");if(!stacked)ImGui::SameLine();
                if(ImGui::Button(fit_label)){
                    double norm2=1;
                    if(inspection && !inspection->selection.normalize && !inspection->metric_norm2.empty())
                        norm2=*std::max_element(inspection->metric_norm2.begin(),inspection->metric_norm2.end());
                    isovalue=std::clamp(static_cast<float>(.03*std::sqrt(std::max(0.0,norm2))),1e-7f,.2f);
                }
                cov::validation::item("scene.fit_component");
                ImGui::PushStyleColor(ImGuiCol_Text,ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
                ImGui::TextWrapped("%s",threshold_help);
                ImGui::PopStyleColor();
                ImGui::Checkbox(cov::ui::aomo_text(language,"Show numbers"),&molecule_render.show_numbers);
                cov::validation::item("scene.numbers");
                if(molecule_render.show_numbers) {
                    ImGui::Checkbox(cov::ui::aomo_text(language,"Atom numbers"),&molecule_render.number_atoms);
                    cov::validation::item("scene.numbers.atoms");
                    ImGui::Checkbox(cov::ui::aomo_text(language,"Ligand numbers"),&molecule_render.number_fragments);
                    cov::validation::item("scene.numbers.fragments");
                    ImGui::Checkbox(cov::ui::aomo_text(language,"Skip H labels"),&molecule_render.number_ignore_h);
                    cov::validation::item("scene.numbers.ignore_h");
                }
                ImGui::End();
            }

            ImGui::SetNextWindowPos(ImVec2(layout.controls.x, layout.controls.y), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(layout.controls.width, layout.controls.height), ImGuiCond_Always);
            ImGui::SetNextWindowBgAlpha(0.965f);
            constexpr ImGuiWindowFlags panel_flags =
                ImGuiWindowFlags_NoTitleBar |
                ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoSavedSettings |
                cov::ui::background_panel_flags;

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
                cov::validation::field("status.detail",status_detail);
            } else {
                disabled_wrapped(cov::ui::tr(cov::ui::Text::IdleHint, language));
            }
            if(status==StatusKind::Error && !status_error_detail.empty() &&
               ImGui::TreeNode(scene_text(language,"Error details","错误详情","エラーの詳細","Détails de l’erreur"))) {
                disabled_wrapped(status_error_detail.c_str());
                ImGui::TreePop();
            }
            ImGui::Separator();
            ImGui::Spacing();

            const bool show_input_panels=!aomo_available || ImGui::CollapsingHeader(
                scene_text(language,"Input files and calculation details","输入文件与计算详情",
                    "入力ファイルと計算の詳細","Fichiers et détails du calcul"));
            cov::validation::item("input.details");
            if(show_input_panels){
            cov::ui::begin_card("##file_card", 192.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::FileSection, language));
            std::optional<std::filesystem::path> recent_to_load;
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::OpenFile, language),
                              ImVec2(150.0f * ui_scale, 0.0f))) {
                const cov::FileDialogResult dialog = cov::open_wavefunction_file_dialog(language);
                if (dialog.selected()) {
                    load_inputs({dialog.path});
                } else if (!dialog.supported) {
                    ImGui::OpenPopup("##file_browser");
                } else if (!dialog.cancelled && !dialog.error.empty()) {
                    status = StatusKind::Error;
                    status_detail = dialog.supported
                                        ? dialog.error
                                        : cov::ui::tr(cov::ui::Text::OpenDialogUnsupported, language);
                }
            }
            if (const auto selected = file_browser(language)) load_inputs({*selected});
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
                load_inputs({path_from_utf8(path_buffer.data())});
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
            if (recent_to_load) load_inputs({*recent_to_load});
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
            cov::validation::item("input.frame_tracking");
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
                    if(frame_tracking->tracking_budget_exhausted){
                        const std::string unresolved=std::to_string(frame_tracking->unresolved_from.size())+
                            " / "+std::to_string(frame_tracking->unresolved_to.size());
                        metric_row(cov::ui::aomo_text(language,"Unresolved correspondence"),unresolved.c_str());
                        metric_row(cov::ui::tr(cov::ui::Text::TrackingOptimisation,language),
                            cov::ui::aomo_text(language,"Matching limit reached"));
                    }else{
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
                    }
                    ImGui::EndTable();
                }
            } else {
                ImGui::TextDisabled("%s", cov::ui::tr(
                    cov::ui::Text::NoPreviousFrame, language));
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            }
            {
                orbital_ui.active_view=active_view();
                nbo_ui.aomo.active_view=orbital_ui.active_view;
                orbital_ui.inspection=inspection.get();
                cov::ui::OrbitalUIActions orbital_actions;
                const bool show_browser=!aomo_available || ImGui::CollapsingHeader(
                    scene_text(language,"Full MO browser","完整 MO 浏览器","全 MO 一覧","Liste complète des OM"));
                cov::validation::item("browser.expand");
                if(show_browser){
                cov::validation::anchor("panel.browser");
                cov::ui::begin_card("##orbital_browser_card", 0);
                cov::ui::section_title(cov::ui::tr(cov::ui::Text::OrbitalBrowser, language));
                if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                    cov::ui::draw_orbital_browser(*wavefunction, canonical_mo_index, orbital_ui,
                                                  language, ui_scale, orbital_actions);
                } else {
                    ImGui::TextDisabled("—");
                }
                cov::ui::end_card();
                ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

                }

                cov::validation::anchor("panel.diagram");
                cov::ui::begin_card("##energy_diagram_card", 0);
                cov::ui::section_title(cov::ui::tr(cov::ui::Text::EnergyDiagram, language));
                cov::ui::OrbitalUIActions diagram_actions;
                if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                    if(integration)cov::ui::prepare_nbo_aomo_state(nbo_ui.aomo,*integration,*wavefunction);
                    const bool atom_numbers=molecule_render.show_numbers && molecule_render.number_atoms;
                    const bool fragment_numbers=molecule_render.show_numbers && molecule_render.number_fragments;
                    if(nbo_ui.aomo.show_atom_numbers!=atom_numbers ||
                       nbo_ui.aomo.show_fragment_numbers!=fragment_numbers ||
                       nbo_ui.aomo.number_ignore_h!=molecule_render.number_ignore_h)++nbo_ui.aomo.revision;
                    nbo_ui.aomo.show_atom_numbers=atom_numbers;
                    nbo_ui.aomo.show_fragment_numbers=fragment_numbers;
                    nbo_ui.aomo.number_ignore_h=molecule_render.number_ignore_h;
                    cov::ui::draw_energy_diagram(*wavefunction, canonical_mo_index, orbital_ui,
                                                 language, ui_scale, diagram_actions);
                } else {
                    ImGui::TextDisabled("—");
                }
                cov::ui::end_card();
                ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

                if (orbital_actions.select_orbital) pending_mo_index = orbital_actions.select_orbital;
                if (diagram_actions.select_orbital) pending_mo_index = diagram_actions.select_orbital;
                canonical_requested=orbital_actions.select_orbital.has_value() || diagram_actions.select_orbital.has_value();
                if (nbo_ui.focus.pending_canonical_selection) {
                    pending_mo_index = nbo_ui.focus.pending_canonical_selection;
                    canonical_requested=true;
                    nbo_ui.focus.pending_canonical_selection.reset();
                }
                const bool export_analysis=orbital_actions.export_analysis || diagram_actions.export_analysis;
                const bool export_requested=orbital_actions.export_diagram || diagram_actions.export_diagram || export_analysis;
                const auto export_content=export_analysis?cov::DiagramExportContent::AnalysisData:cov::DiagramExportContent::Images;
                if (export_requested && wavefunction) {
                    std::filesystem::path base = current_file.empty()
                                                     ? std::filesystem::current_path() / "mo_diagram"
                                                     : current_file;
                    base = cov::validation::export_base(base);
                    const auto snapshot=diagram_actions.drawn_diagram;
                    cov::MODiagramExportResult result;
                    if (snapshot) {
                        auto presentation=snapshot->options;
                        presentation.display_names.clear();
                        presentation.display_irreps.clear();
                        presentation.display_point_groups.clear();
                        presentation.display_names.reserve(wavefunction->orbitals.size());
                        presentation.display_irreps.reserve(wavefunction->orbitals.size());
                        presentation.display_point_groups.reserve(wavefunction->orbitals.size());
                        const auto names=nbo_ui.aomo.names?nbo_ui.aomo.names:cov::ui::canonical_mo_names(*wavefunction);
                        for(std::size_t i=0;i<wavefunction->orbitals.size();++i) {
                            const auto* name=names&&i<names->canonical.size()?&names->canonical[i]:nullptr;
                            presentation.display_names.push_back(cov::ui::canonical_mo_display_label(*wavefunction,i,name));
                            presentation.display_irreps.push_back(cov::ui::canonical_mo_current_irrep(*wavefunction,i,name));
                            presentation.display_point_groups.push_back(name&&name->verified?name->point_group:std::string{});
                        }
                        presentation.figure_title=scene_text(language,"Molecular orbital energies","分子轨道能级",
                            "分子軌道のエネルギー","Énergies des orbitales moléculaires");
                        presentation.axis_title=presentation.energy_axis_mode==cov::EnergyAxisMode::Linear
                            ? scene_text(language,"Energy","能量","エネルギー","Énergie")
                            : scene_text(language,"Nonlinear energy axis","非线性能量轴",
                                "非線形エネルギー軸","Axe d’énergie non linéaire");
                        presentation.raster_text=cov::ui::raster_text;
                        presentation.raster_text_width=cov::ui::raster_text_width;
                        result=cov::export_mo_diagram_bundle({snapshot->data,std::move(presentation)},base,export_content);
                    }
                    else result.error="No current diagram view is available for export";
                    const bool exported=export_analysis?(result.json&&result.csv):(result.svg&&result.png);
    #ifdef COV_ENABLE_VALIDATION
                    cov::validation::record("export.actual","{\"base\":"+cov::validation::quote(path_to_utf8(base))+
                        ",\"snapshot_id\":"+(snapshot?cov::validation::quote(snapshot->data.view->id):"null")+
                        ",\"mode\":"+(snapshot?std::to_string(static_cast<int>(snapshot->data.mode)):"null")+
                        ",\"selected_index\":"+(snapshot && snapshot->data.view->inspected_orbital_index
                            ?std::to_string(*snapshot->data.view->inspected_orbital_index):"null")+
                        ",\"content\":"+cov::validation::quote(export_analysis?"analysis_data":"images")+
                        ",\"success\":"+(exported?"true":"false")+"}");
    #endif
                    if (exported) {
                        if(export_analysis) export_analysis_companions(base);
                        status = StatusKind::Exported;
                        status_detail = path_to_utf8(result.svg_path.parent_path() /
                            result.svg_path.stem()) + (export_analysis?".{json,csv}":".{png,svg}");
                    } else {
                        status = StatusKind::Error;
                        status_detail=cov::ui::tr(cov::ui::Text::ExportFailed,language);
                        status_error_detail=result.error;
                    }
                }
            }

            cov::validation::anchor("panel.nbo");
            cov::ui::begin_card("##nbo_card", 0);
            const auto nbo_actions = cov::ui::draw_nbo_panel(
                nbo_ui, language, static_cast<bool>(wavefunction),
                nbo_wavefunction.has_value(), active_view().kind, ui_scale,
                wavefunction ? &*wavefunction : nullptr, canonical_mo_index,
                orbital_ui.diagram_cache.snapshot.get());
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));
            const bool attach_requested = nbo_actions.attach;
            if(nbo_actions.choose_input)choose_nbo_input=true;
            std::optional<bool> requested_set;
            if (nbo_actions.canonical_set) requested_set=false;
            if (nbo_actions.nbo_set) requested_set=true;
            if(requested_set && *requested_set && nbo_mo_index==std::numeric_limits<std::size_t>::max()){
                requested_set.reset();nbo_ui.error="Select a specific available NBO orbital first";
            }
            if (nbo_actions.selected_orbital && nbo_ui.dataset) {
                nbo_ui.selected_orbital=*nbo_actions.selected_orbital;
                nbo_mo_index=nbo_ui.selected_orbital;
                if (nbo_active) pending_mo_index=nbo_mo_index;
            }
            if (nbo_actions.export_bundle && nbo_ui.dataset) {
                try {
                    auto base=nbo_ui.export_path[0]
                        ? path_from_utf8(nbo_ui.export_path.data())
                        : path_from_utf8(nbo_ui.dataset->source.path);
                    base=cov::validation::export_base(base);
                    cov::ui::export_nbo_bundle(*nbo_ui.dataset,nbo_ui.selected_orbital,base,
                                               wavefunction ? &*wavefunction : nullptr,
                                               canonical_mo_index,nbo_ui.contribution_threshold,
                                               active_view(),integration?&*integration:nullptr,
                                               orbital_ui.diagram_cache.snapshot.get(),&nbo_ui.focus,language);
                    export_analysis_companions(base);
                    nbo_ui.export_status=path_to_utf8(base)+".{nbo.json,npa.csv,nao.csv,nbo.csv,wiberg.csv,e2.csv,e2-sections.csv,view.json,view.svg,view.png,focus.json,focus.csv,focus.groups.csv,focus.svg,focus.png}";
                    cov::validation::record("nbo.export","{\"base\":"+
                        cov::validation::quote(path_to_utf8(base))+
                        ",\"selected_index\":"+std::to_string(nbo_ui.selected_orbital)+
                        ",\"source\":"+cov::validation::quote(nbo_ui.dataset->source.path)+"}");
                } catch (const std::exception& e) { nbo_ui.export_status=e.what(); }
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
                if (fields_ready()) cov::ui::status_badge(cov::ui::tr(cov::ui::Text::GPUResident, language), cov::ui::Tone::Success);
                else if (fields_busy()) disabled_wrapped(status_label(StatusKind::Computing, language));
            } else {
                ImGui::TextDisabled("%s —", cov::ui::tr(cov::ui::Text::CUDADevice, language));
            }
            ImGui::TextDisabled("%s: %s",
                                cov::ui::tr(cov::ui::Text::FontStatus, language),
                                cov::ui::font_status(language));
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::InteractionHint, language));
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::IsovalueHint, language));
            cov::ui::end_card();

            ImGui::EndChild();
            ImGui::End();
            cov::validation::field("language",std::to_string(static_cast<int>(language)));
            identity();
            cov::validation::ui_frame(mo_index,pending_mo_index.value_or(mo_index));

            if(cov::validation::forensic_mode())
                cov::validation::record("forensic.input","{\"canonical_path\":"+
                    cov::validation::quote(path_to_utf8(current_file))+
                    ",\"status\":"+std::to_string(static_cast<int>(status))+
                    ",\"error\":"+cov::validation::quote(status_error_detail)+
                    ",\"nbo_source\":"+cov::validation::quote(integration?integration->dataset.source.path:std::string{})+"}");
            ImGui::Render();
            ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
            cov::validation::end_frame(fb_w,fb_h,mo_index,orbital_ui,active_wavefunction());
            if(profile_open && profile_draw_frames<3) {
                unsigned char pixel[4]{};
                glReadPixels(20,std::max(0,fb_h-20),1,1,GL_RGBA,GL_UNSIGNED_BYTE,pixel);
                std::fprintf(stderr,"COV initial draw: framebuffer=%dx%d vertices=%d rgba=%u,%u,%u,%u gl_error=%u renderer=%s\n",
                    fb_w,fb_h,ImGui::GetDrawData()->TotalVtxCount,pixel[0],pixel[1],pixel[2],pixel[3],
                    glGetError(),reinterpret_cast<const char*>(glGetString(GL_RENDERER)));
            }
            glfwSwapBuffers(window);
            if(profile_open && profile_draw_frames<3) {
                GLint saved_read=0;glGetIntegerv(GL_READ_BUFFER,&saved_read);glReadBuffer(GL_FRONT);
                unsigned char pixel[4]{};glReadPixels(50,std::max(0,fb_h-50),1,1,GL_RGBA,GL_UNSIGNED_BYTE,pixel);
                glReadBuffer(saved_read);
                std::fprintf(stderr,"COV initial swap: frame=%u front=%u,%u,%u,%u\n",profile_draw_frames,pixel[0],pixel[1],pixel[2],pixel[3]);
#ifdef _WIN32
                std::fprintf(stderr,"COV present window: hwnd=%p dc=%p\n",WindowFromDC(wglGetCurrentDC()),wglGetCurrentDC());
#endif
                ++profile_draw_frames;
            }
            if(profile_frame_pending){profile_stage("first-frame");profile_frame_pending=false;}
            if (attach_requested) {
                try {
                    cov::NboReadOptions options;
                    if (nbo_ui.analysis_segment >= 0)
                        options.analysis_segment = static_cast<std::size_t>(nbo_ui.analysis_segment);
                    if (nbo_ui.archive47[0]) options.archive47=path_from_utf8(nbo_ui.archive47.data());
                    if (nbo_ui.aonbo[0]) options.aonbo=path_from_utf8(nbo_ui.aonbo.data());
                    if (nbo_ui.nbomo[0]) options.nbomo=path_from_utf8(nbo_ui.nbomo.data());
                    if (nbo_ui.naomo[0]) options.naomo=path_from_utf8(nbo_ui.naomo.data());
                    if (nbo_ui.aonao[0]) options.aonao=path_from_utf8(nbo_ui.aonao.data());
                    if (nbo_ui.naonbo[0]) options.naonbo=path_from_utf8(nbo_ui.naonbo.data());
                    auto dataset=cov::read_nbo(path_from_utf8(nbo_ui.path.data()),options);
                    attach_integration(cov::integrate_nbo(*wavefunction,dataset));
                } catch (const std::exception& e) {
                    nbo_ui.error=e.what();
                    cov::validation::record("nbo.attach.error","{\"reason\":"+
                        cov::validation::quote(nbo_ui.error)+"}");
                }
            }

            if(nbo_ui.pending_candidate && nbo_ui.input_discovery){
                const auto selected=*nbo_ui.pending_candidate;nbo_ui.pending_candidate.reset();
                if(selected<nbo_ui.input_discovery->candidates.size()){
                    const auto candidate=nbo_ui.input_discovery->candidates[selected];
                    try {apply_candidate(candidate);}catch(const std::exception& e){
                        clear_integration();nbo_ui.input_discovery.reset();nbo_ui.error=e.what();
                        nbo_ui.input_status="NBO association failed; Gaussian view is retained";
                        status=StatusKind::Error;status_detail=nbo_ui.error;
                    }
                }
            }
            if(nbo_ui.aomo.pending_selection){
                const auto selection=*nbo_ui.aomo.pending_selection;nbo_ui.aomo.pending_selection.reset();
                try{apply_selection(selection);}catch(const std::exception& e){nbo_ui.aomo.status=e.what();
                    nbo_ui.error=e.what();status=StatusKind::Error;status_detail=e.what();
                    cov::validation::record("aomo.selection.error","{\"reason\":"+cov::validation::quote(e.what())+"}");}
            }
            if(nbo_ui.aomo.export_requested){
                nbo_ui.aomo.export_requested=false;
                if(integration && nbo_ui.aomo.drawn_snapshot){
                    const auto base=cov::validation::export_base(nbo_ui.aomo.export_path[0]?
                        path_from_utf8(nbo_ui.aomo.export_path.data()):current_file);
                    const auto content=nbo_ui.aomo.export_content;
                    const auto result=cov::ui::export_nbo_aomo_bundle(*nbo_ui.aomo.drawn_snapshot,*integration,base,content);
                    if(content==cov::DiagramExportContent::AnalysisData && result.json && result.csv)
                        export_analysis_companions(base);
                    nbo_ui.aomo.export_status=result.error.empty()?path_to_utf8(
                        content==cov::DiagramExportContent::AnalysisData?result.json_path:result.svg_path):result.error;
                }
            }
            if(canonical_requested){
                const auto selected=pending_mo_index;
                try{activate_set(false);pending_mo_index=selected;}
                catch(const std::exception& e){nbo_ui.error=e.what();}
            }
            if (requested_set) {
                try { activate_set(*requested_set); }
                catch (const std::exception& e) { nbo_ui.error=e.what(); }
            }
            // Selection debounce: at most the latest requested orbital is evaluated
            // once at the end of this frame. Browser hover/filtering does not recompute the grid.
            if (pending_mo_index && active_wavefunction() &&
                *pending_mo_index < active_wavefunction()->orbitals.size()) {
                if (*pending_mo_index != mo_index) {
                    mo_index = *pending_mo_index;
                    if (nbo_active) nbo_mo_index=mo_index;
                    else canonical_mo_index=mo_index;
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

            if (cov::validation::done()) {exit_code=cov::validation::result();break;}
        }

        if (evaluator) evaluator->detach_gl_texture();
    } catch (const std::exception& e) {
        report_fatal_error(e.what());
        exit_code = 1;
    }

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return exit_code;
}
