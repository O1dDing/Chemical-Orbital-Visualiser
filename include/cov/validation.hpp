#pragma once
#include "cov/orbital_ui.hpp"
#include "cov/volume_renderer.hpp"
#include "cov/viewer_layout.hpp"
#include <imgui.h>
#include <filesystem>
#include <string>

namespace cov::validation {
#ifdef COV_ENABLE_VALIDATION
bool configure(int argc, char** argv);
bool active();
bool forensic_mode();
bool background();
int window_width();
int window_height();
bool done();
int result();
void begin_frame(OrbitCamera&, MoleculeRenderSettings&, float&, int&, bool&);
void input_frame();
std::vector<std::filesystem::path> take_dropped_paths();
void evaluated(std::size_t mo, const char* reason, float milliseconds);
void orbital_identity(const std::string& orbital_set, const std::string& dataset,
                      const std::string& spin, std::size_t source_index,
                      const std::string& association,
                      const std::string& coefficient_source,
                      bool direct_fchk_coefficients, bool density_verified);
void ui_frame(std::size_t drawn, std::size_t requested);
void after_scene(const VolumeRenderer&, const GridBox&, std::size_t mo, std::size_t field_index=0);
void scene_view(const ViewerLayout&, const OrbitCamera&);
void end_frame(int width, int height, std::size_t applied,
               const ui::OrbitalUIState&, const Wavefunction*);
void item(const std::string& id);
void hit(const std::string& id, ImVec2 lo, ImVec2 hi);
void chrome_hit(const std::string& id, ImVec2 lo, ImVec2 hi);
void anchor(const std::string& id);
void record(const std::string& kind, const std::string& json);
void field(const std::string& label, const std::string& value);
std::string quote(const std::string& value);
std::filesystem::path export_base(const std::filesystem::path& original);
#else
inline bool configure(int, char**) { return false; }
inline bool active() { return false; }
inline bool forensic_mode() { return false; }
inline bool background() { return false; }
inline int window_width() { return 2100; }
inline int window_height() { return 1250; }
inline bool done() { return false; }
inline int result() { return 0; }
inline void begin_frame(OrbitCamera&, MoleculeRenderSettings&, float&, int&, bool&) {}
inline void input_frame() {}
inline std::vector<std::filesystem::path> take_dropped_paths() {return {};}
inline void evaluated(std::size_t, const char*, float) {}
inline void orbital_identity(const std::string&, const std::string&,
                             const std::string&, std::size_t,
                             const std::string&, const std::string&, bool, bool) {}
inline void ui_frame(std::size_t, std::size_t) {}
inline void after_scene(const VolumeRenderer&, const GridBox&, std::size_t, std::size_t=0) {}
inline void scene_view(const ViewerLayout&, const OrbitCamera&) {}
inline void end_frame(int, int, std::size_t, const ui::OrbitalUIState&, const Wavefunction*) {}
inline void item(const std::string&) {}
inline void hit(const std::string&, ImVec2, ImVec2) {}
inline void chrome_hit(const std::string&, ImVec2, ImVec2) {}
inline void anchor(const std::string&) {}
inline void record(const std::string&, const std::string&) {}
inline void field(const std::string&, const std::string&) {}
inline std::string quote(const std::string&) { return {}; }
inline std::filesystem::path export_base(const std::filesystem::path& p) { return p; }
#endif
}
