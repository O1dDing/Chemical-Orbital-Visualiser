#include "cov/validation.hpp"
#include "cov/validation_navigation.hpp"
#include "cov/forensic_capture.hpp"
#include "cov/gl_api.hpp"
#include <imgui_internal.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <map>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace cov::validation {
namespace {
using Clock = std::chrono::steady_clock;
struct Target { ImVec2 lo, hi; ImGuiWindow* window; ImRect clip; };
struct Command { std::string op, id, value; std::vector<float> args; std::vector<std::string> paths; };
std::vector<std::filesystem::path> dropped_paths;
bool enabled = false;
bool hidden_window = false;
bool forensic = false;
int detail_page = 0;
float detail_scroll_before = -1;
bool detail_page_pending = false;
int requested_window_width = 2100, requested_window_height = 1250;
std::filesystem::path output;
std::string export_name="actual-export";
std::vector<Command> commands;
std::size_t next = 0, frame = 0, generation = 0, rendered_generation = 0;
std::size_t volume_mo = 0, rendered_mo = 0;
std::size_t drawn_ui_mo = 0, requested_mo = 0, diagram_generation = 0;
std::string active_set="canonical", active_dataset, active_spin="alpha", active_association="not_attached";
std::string active_coefficient_source="not_available";
bool active_direct_fchk=false, active_density_verified=false;
std::size_t active_source_index=std::numeric_limits<std::size_t>::max();
std::string rendered_set="canonical", rendered_dataset, rendered_spin="alpha", rendered_association="not_attached";
std::string rendered_coefficient_source="not_available";
bool rendered_direct_fchk=false, rendered_density_verified=false;
std::size_t rendered_source_index=std::numeric_limits<std::size_t>::max();
int stage = 0, attempts = 0, failures = 0, cooldown = 0;
bool complete_command = false;
Clock::time_point started, command_started;
std::ofstream frames, actions, events;
std::map<std::string, Target> targets, previous;
std::vector<std::string> trace;
std::string evaluation_reason;
float kernel_ms = 0;
ImVec2 injected_mouse(-100,-100);
std::string scene_view_json = "null";
bool volume_command(const Command& c) { return c.op=="volume" || c.op=="volume_full"; }
void write_volume_binary(const std::filesystem::path& path, const std::vector<float>& volume) {
    static_assert(sizeof(float)==4 && std::numeric_limits<float>::is_iec559,
                  "validation texture evidence requires IEEE-754 float32");
    std::ofstream file(path,std::ios::binary);
    if(!file)throw std::runtime_error("cannot create complete texture evidence");
    const std::uint32_t probe=1;
    if(*reinterpret_cast<const unsigned char*>(&probe)==1) {
        file.write(reinterpret_cast<const char*>(volume.data()),
                   static_cast<std::streamsize>(volume.size()*sizeof(float)));
    } else {
        std::vector<unsigned char> bytes(volume.size()*4);
        for(std::size_t i=0;i<volume.size();++i) {
            std::uint32_t bits;std::memcpy(&bits,&volume[i],4);
            for(int j=0;j<4;++j)bytes[4*i+j]=static_cast<unsigned char>((bits>>(8*j))&255);
        }
        file.write(reinterpret_cast<const char*>(bytes.data()),static_cast<std::streamsize>(bytes.size()));
    }
    file.close();
    if(!file)throw std::runtime_error("complete texture evidence write failed");
}
double elapsed() { return std::chrono::duration<double>(Clock::now()-started).count(); }
std::string number(double v) { std::ostringstream s; s << std::setprecision(17) << v; return s.str(); }
void finish(const std::string& status, const std::string& detail = {}) {
    const auto& c = commands[next];
    actions << "{\"schema\":1,\"command\":" << next << ",\"op\":" << quote(c.op)
            << ",\"id\":" << quote(c.id) << ",\"status\":" << quote(status)
            << ",\"detail\":" << quote(detail) << ",\"frame\":" << frame
            << ",\"seconds\":" << number(std::chrono::duration<double>(Clock::now()-command_started).count()) << "}\n";
    actions.flush();
    if (status != "executed") ++failures;
    ++next; stage = attempts = cooldown = detail_page = 0; complete_command = false;
    detail_scroll_before=-1;detail_page_pending=false;
    command_started = Clock::now();
}
NavigationTarget navigation_target(const Target& t) { return {t.lo,t.hi,t.window,t.clip.Min,t.clip.Max}; }
bool point_visible(const Target& t) {
    return navigation_target_visible(navigation_target(t));
}
bool seek(const Target& t, bool reveal_entire_item = false) {
    auto& io=ImGui::GetIO();
    const auto step=plan_navigation(navigation_target(t),reveal_entire_item);
    if (step.kind==NavigationKind::Unreachable) return false;
    injected_mouse=step.mouse;io.AddMousePosEvent(step.mouse.x,step.mouse.y);
    if (step.kind==NavigationKind::Ready) return true;
    events<<"{\"kind\":\"input.seek\",\"frame\":"<<frame<<",\"command\":"<<next
          <<",\"phase\":"<<quote(step.kind==NavigationKind::Move?"move":"wheel")
          <<",\"window\":"<<quote(step.scrolling_window?step.scrolling_window->Name:"")
          <<",\"mouse\":["<<step.mouse.x<<','<<step.mouse.y<<"],\"wheel\":["
          <<step.wheel.x<<','<<step.wheel.y<<"],\"target\":["<<t.lo.x<<','<<t.lo.y<<','<<t.hi.x<<','<<t.hi.y
          <<"],\"clip\":["<<t.clip.Min.x<<','<<t.clip.Min.y<<','<<t.clip.Max.x<<','<<t.clip.Max.y<<"]}\n";
    if (step.kind==NavigationKind::Wheel) {
        io.AddMouseWheelEvent(step.wheel.x,step.wheel.y);cooldown=3;
    }
    return false;
}
void framebuffer(const std::filesystem::path& path, int w, int h) {
    // Read the completed, visible, double-buffered viewer + ImGui client area.
    GLint read_buffer=0, alignment=0;
    glGetIntegerv(GL_READ_BUFFER,&read_buffer); glGetIntegerv(GL_PACK_ALIGNMENT,&alignment);
    glReadBuffer(GL_BACK); glPixelStorei(GL_PACK_ALIGNMENT,1);
    std::vector<unsigned char> rgba(static_cast<std::size_t>(w)*h*4);
    glReadPixels(0,0,w,h,GL_RGBA,GL_UNSIGNED_BYTE,rgba.data());
    glReadBuffer(read_buffer); glPixelStorei(GL_PACK_ALIGNMENT,alignment);
    if (glGetError()!=GL_NO_ERROR) throw std::runtime_error("framebuffer readback failed");
    for (std::size_t i=0;i<rgba.size();i+=4) std::swap(rgba[i],rgba[i+2]);
    std::ofstream s(path,std::ios::binary);
    auto u16=[&](std::uint16_t v){s.write(reinterpret_cast<char*>(&v),2);};
    auto u32=[&](std::uint32_t v){s.write(reinterpret_cast<char*>(&v),4);};
    s.write("BM",2); u32(54+static_cast<std::uint32_t>(rgba.size())); u32(0); u32(54);
    u32(40); u32(w); u32(h); u16(1); u16(32); u32(0); u32(static_cast<std::uint32_t>(rgba.size()));
    u32(0);u32(0);u32(0);u32(0); s.write(reinterpret_cast<char*>(rgba.data()),rgba.size());
    if (!s) throw std::runtime_error("cannot write framebuffer");
}
std::string state_json(std::size_t applied, const ui::OrbitalUIState& ui, const Wavefunction* wf) {
    std::ostringstream s; s << std::setprecision(17);
    s << "{\"schema\":1,\"frame\":" << frame << ",\"elapsed_seconds\":" << elapsed()
      << ",\"rendered_mo\":" << rendered_mo << ",\"applied_mo\":" << applied
      << ",\"rendered_set\":" << quote(rendered_set)
      << ",\"rendered_dataset\":" << quote(rendered_dataset)
      << ",\"rendered_spin\":" << quote(rendered_spin)
      << ",\"rendered_source_index\":";
    if(rendered_source_index==std::numeric_limits<std::size_t>::max())s<<"null";else s<<rendered_source_index;
    s << ",\"applied_set\":" << quote(active_set)
      << ",\"applied_dataset\":" << quote(active_dataset)
      << ",\"applied_spin\":" << quote(active_spin)
      << ",\"applied_source_index\":";
    if(active_source_index==std::numeric_limits<std::size_t>::max())s<<"null";else s<<active_source_index;
    s << ",\"association\":" << quote(active_association)
      << ",\"coefficient_source\":" << quote(active_coefficient_source)
      << ",\"direct_fchk_coefficients\":" << (active_direct_fchk?"true":"false")
      << ",\"density_verified\":" << (active_density_verified?"true":"false")
      << ",\"rendered_coefficient_source\":" << quote(rendered_coefficient_source)
      << ",\"rendered_direct_fchk_coefficients\":" << (rendered_direct_fchk?"true":"false")
      << ",\"rendered_density_verified\":" << (rendered_density_verified?"true":"false")
      << ",\"drawn_ui_mo\":" << drawn_ui_mo << ",\"requested_mo\":" << requested_mo
      << ",\"scene_view\":" << scene_view_json
      << ",\"diagram_generation\":" << diagram_generation
      << ",\"rendered_generation\":" << rendered_generation << ",\"volume_generation\":" << generation
      << ",\"scene_matches_applied\":" << (rendered_mo==applied && rendered_set==active_set &&
          rendered_dataset==active_dataset && rendered_spin==active_spin &&
          rendered_source_index==active_source_index && rendered_coefficient_source==active_coefficient_source?"true":"false")
      << ",\"evaluation_reason\":" << quote(evaluation_reason) << ",\"kernel_ms\":" << kernel_ms
      << ",\"compact\":" << (ui.hide_ligand_centred_intermediates?"true":"false")
      << ",\"energy_unit\":" << static_cast<int>(ui.energy_unit)
      << ",\"axis_mode\":" << static_cast<int>(ui.energy_axis_mode)
      << ",\"filter\":" << static_cast<int>(ui.filter.mode);
    if (wf && applied<wf->orbitals.size()) {
        const auto& mo=wf->orbitals[applied];
        s << ",\"energy_hartree\":";
        if(std::isfinite(mo.energy_hartree) && active_set=="canonical")s<<mo.energy_hartree;
        else s<<"null";
        s << ",\"energy_semantics\":" << quote(active_set=="canonical"?"canonical eigenvalue":"not applicable; NBO diagonal Fock is separate")
          << ",\"occupation\":";
        if(active_set=="canonical")s<<mo.occupation;else s<<"null";
        s << ",\"occupation_semantics\":" << quote(active_set=="canonical"?"canonical occupation":"selected orbital or combination metadata; renderer placeholder is not physical occupation")
          << ",\"spin\":";
        if(active_set=="canonical")s<<static_cast<int>(mo.spin);else s<<"null";
        s << ",\"spin_semantics\":" << quote(active_set=="canonical"?"canonical wavefunction spin":"NBO spin is applied_spin; renderer channel may differ");
    }
    s << '}'; return s.str();
}
}

std::string quote(const std::string& v) {
    std::string s="\"";
    for (unsigned char c:v) {
        if(c=='"'||c=='\\') {s+='\\';s+=c;}
        else if(c=='\n') s+="\\n";
        else if(c=='\r') s+="\\r";
        else if(c=='\t') s+="\\t";
        else if(c<32) {const char* h="0123456789abcdef";s+="\\u00";s+=h[c>>4];s+=h[c&15];}
        else s+=c;
    }
    return s+'"';
}
bool configure(int argc, char** argv) {
    std::filesystem::path plan;
    bool input_seen = false;
    for(int i=1;i<argc;++i) {
        const std::string a=argv[i];
        if(a.starts_with("--compute-backend=") || a.starts_with("--compute-device=")) continue;
        if(!a.starts_with("--") && !input_seen) { input_seen=true; continue; }
        if(a=="--validation-plan" && i+1<argc) plan=std::filesystem::u8path(argv[++i]);
        else if(a=="--validation-output" && i+1<argc) output=std::filesystem::u8path(argv[++i]);
        else if(a=="--validation-background") hidden_window=true;
        else if(a=="--validation-forensic") forensic=true;
        else throw std::runtime_error("unknown/incomplete validation argument: "+a);
    }
    if(plan.empty() && output.empty() && !hidden_window && !forensic) return false;
    if(plan.empty() || output.empty()) throw std::runtime_error("plan and output are both required");
    std::ifstream in(plan); std::string line;
    std::getline(in,line); if(line!="COV_VALIDATION 1") throw std::runtime_error("unsupported validation plan schema");
    while(std::getline(in,line)) {
        if(line.empty() || line[0]=='#') continue;
        std::istringstream r(line); Command c; r>>c.op;
        if(c.op=="drop") {std::string path;while(r>>std::quoted(path))c.paths.push_back(path);if(c.paths.empty())throw std::runtime_error("drop requires one or more quoted paths");}
        else if(c.op=="scene") {float x;while(r>>x)c.args.push_back(x);if(c.args.size()!=6)throw std::runtime_error("scene requires opacity yaw pitch distance iso resolution");}
        else if(c.op=="window") {
            float width=0,height=0; std::string extra;
            if(!(r>>width>>height) || (r>>extra) || !std::isfinite(width) || !std::isfinite(height) ||
               width<640 || width>7680 || height<360 || height>4320 ||
               std::floor(width)!=width || std::floor(height)!=height)
                throw std::runtime_error("window requires integral width 640..7680 and height 360..4320");
            c.args={width,height};
        }
        else {
            r>>std::quoted(c.id);
            if(c.op=="text" || volume_command(c)) r>>std::quoted(c.value);
            if(c.op=="drag" || c.op=="wheel") {
                float x=0,y=0; std::string extra;
                if(!(r>>x) || (c.op=="drag" && !(r>>y)) || (r>>extra) ||
                   !std::isfinite(x) || !std::isfinite(y) || std::abs(x)>10000 || std::abs(y)>10000)
                    throw std::runtime_error("drag/wheel requires finite bounded pointer deltas");
                c.args={x,y};
            }
        }
        if(c.op!="drop"&&c.op!="scene"&&c.op!="window"&&c.op!="drag"&&c.op!="wheel"&&c.op!="click"&&c.op!="expand"&&c.op!="hover"&&c.op!="seek"&&c.op!="text"&&c.op!="capture"&&!volume_command(c)&&c.op!="key"&&c.op!="wait"&&c.op!="export-name"&&c.op!="inspect"&&c.op!="inspect-details") throw std::runtime_error("unknown plan command");
        if((c.op=="inspect" || c.op=="inspect-details") && !forensic)
            throw std::runtime_error("inspect requires --validation-forensic");
        if((c.op=="inspect" || c.op=="inspect-details" || c.op=="capture") &&
           (c.id.empty() || c.id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_")!=std::string::npos))
            throw std::runtime_error("capture/inspect requires a plain artifact name");
        if(c.op=="export-name" && (c.id.empty() || c.id.find_first_not_of(
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_")!=std::string::npos)) {
            throw std::runtime_error("export-name requires a plain artifact name");
        }
        commands.push_back(c);
    }
    if(std::filesystem::exists(output / "actions.jsonl")) throw std::runtime_error("refusing to overwrite an existing validation run");
    std::filesystem::create_directories(output);
    std::filesystem::copy_file(plan,output/"plan.txt");
    if(!forensic)frames.open(output/"frames.jsonl");
    actions.open(output/"actions.jsonl");
    events.open(output/"events.jsonl");
    started=command_started=Clock::now(); enabled=true;
    std::ofstream identity(output/"identity.json");
    identity << "{\"schema\":1,\"git_commit\":" << quote(COV_VALIDATION_COMMIT)
             << ",\"input\":" << quote(argv[1]) << ",\"build\":\"validation ON\",\"imgui\":" << quote(IMGUI_VERSION)
             << ",\"window_mode\":" << quote(hidden_window?"background-hidden":"visible")
             << ",\"forensic\":" << (forensic?"true":"false")
             << ",\"protocol\":\"local plan v1\",\"scientific_verdict\":\"external checker required\"}";
    return true;
}
bool active(){return enabled;}
bool forensic_mode(){return enabled&&forensic;}
bool background(){return enabled&&hidden_window;}
int window_width(){return requested_window_width;}
int window_height(){return requested_window_height;}
bool done(){return enabled&&next>=commands.size();}
int result(){return failures?2:0;}
void begin_frame(OrbitCamera& camera, MoleculeRenderSettings& settings, float& iso, int& resolution, bool& resize) {
    if(!enabled)return;
    ++frame; previous=std::move(targets);targets.clear();trace.clear();evaluation_reason.clear();
    if(done())return;
    auto& c=commands[next];
    if(c.op=="window") {
        requested_window_width=static_cast<int>(c.args[0]);
        requested_window_height=static_cast<int>(c.args[1]);
    }
    if(c.op=="scene") {
        settings.orbital_opacity=c.args[0];camera.yaw=c.args[1];camera.pitch=c.args[2];camera.distance=c.args[3];iso=c.args[4];
        const int requested=static_cast<int>(c.args[5]);
        if(requested!=resolution) {resolution=requested;resize=true;}
        complete_command=true;
    }
}
std::vector<std::filesystem::path> take_dropped_paths(){
    auto out=std::move(dropped_paths);dropped_paths.clear();return out;
}
void input_frame() {
    if(!enabled||done())return;
    auto& c=commands[next]; auto& io=ImGui::GetIO();
    // The local plan owns input while active. Discard OS cursor polling queued
    // by the GLFW backend, then feed events through ImGui's normal input API.
    // No Button/Selectable return value or application selection is forced.
    io.ClearEventsQueue();
    io.AddFocusEvent(true);
    io.AddMousePosEvent(injected_mouse.x,injected_mouse.y);
    if(c.op=="drop") {
        if(stage==0)for(const auto& p:c.paths)dropped_paths.push_back(std::filesystem::u8path(p));
        if(++stage>=5)complete_command=true;
        return;
    }
    if(c.op=="scene")return;
    if(c.op=="window") {if(++stage>=4)complete_command=true;return;}
    // Destination naming alone; the following real button click still owns
    // the production export. Existing COV_VALIDATION 1 plans keep the default.
    if(c.op=="export-name") {export_name=c.id;complete_command=true;return;}
    if(volume_command(c)) {++stage;return;}
    if(c.op=="capture" || c.op=="wait" || c.op=="inspect") { if(++stage>=4) complete_command=true;return; }
    if(c.op=="inspect-details") {
        const auto found=previous.find("diagram.details.window");
        if(found==previous.end() || !found->second.window || found->second.window->Collapsed) {
            if(++attempts>=12)finish("failed","details window is not open");
            return;
        }
        auto* w=found->second.window;
        // Real wheel input, never SetScrollY: this also reveals occlusion failures.
        const auto p=ImVec2(w->InnerRect.GetCenter().x,w->InnerRect.GetCenter().y);
        injected_mouse=p;io.AddMousePosEvent(p.x,p.y);
        if(cooldown>0) {io.AddMouseWheelEvent(0,-float(cooldown));cooldown=0;stage=1;return;}
        if(stage==0 && detail_page==0 && w->Scroll.y>1) {
            io.AddMouseWheelEvent(0,100);stage=1;detail_scroll_before=w->Scroll.y;return;
        }
        if(++stage<5)return;
        if(detail_page==0 && w->Scroll.y>1 && detail_scroll_before>=0) {
            finish("failed","details top is unreachable or covered");return;
        }
        if(detail_page>0 && w->Scroll.y<=detail_scroll_before+0.5f && w->Scroll.y<w->ScrollMax.y-1) {
            finish("failed","details scrolling made no progress; check window stacking");return;
        }
        detail_page_pending=true;
        return;
    }
    if(c.op=="key") {
        const std::map<std::string,ImGuiKey> keys={{"Home",ImGuiKey_Home},{"Down",ImGuiKey_DownArrow},{"Up",ImGuiKey_UpArrow},{"Enter",ImGuiKey_Enter},{"Escape",ImGuiKey_Escape}};
        const auto k=keys.find(c.id);if(k==keys.end()){finish("failed","unsupported key");return;}
        if(stage<2)io.AddKeyEvent(k->second,stage==0);
        if(++stage>=4)complete_command=true;return;
    }
    if(cooldown>0){--cooldown;return;}
    if(c.op=="expand" && stage==0) {
        // Observe the registered production header state; never force it.
        // When closed, continue through the ordinary seek/down/up path.
        if(previous.contains(c.id+".open")){finish("executed","header already open");return;}
        if(!previous.contains(c.id)) {
            if(++attempts>=4)finish("executed","header not applicable to this selection");
            return;
        }
        if(!previous.contains(c.id+".closed")){finish("failed","header expansion state not registered");return;}
    }
    if(stage==0) {
        const auto it=previous.find(c.id);
        if(it==previous.end()) {
            if(c.id.rfind("browser.mo.",0)==0) {
                const auto table=previous.find("browser.table");
                if(table!=previous.end()) seek(table->second);
            }
            if(++attempts>60)finish("failed","semantic target not drawn");return;
        }
        if(!seek(it->second,c.op=="hover" || c.op=="seek")) {if(++attempts>60)finish("failed","target clipped or unreachable by wheel input");return;}
        if(c.op=="seek"){complete_command=true;return;}
    }
    // Once the real pointer sequence starts, finish its release and settling
    // frames even if the action closes its own window. Final captures, rather
    // than continued target existence, establish whether the action succeeded.
    if(c.op=="hover") {if(++stage>=4)complete_command=true;return;}
    if(c.op=="wheel" || c.op=="drag") {
        if(c.op=="wheel" && stage==1) io.AddMouseWheelEvent(0,c.args[0]);
        if(c.op=="drag") {
            if(stage==1)io.AddMouseButtonEvent(0,true);
            if(stage==2) {
                injected_mouse.x+=c.args[0];injected_mouse.y+=c.args[1];
                io.AddMousePosEvent(injected_mouse.x,injected_mouse.y);
            }
            if(stage==3)io.AddMouseButtonEvent(0,false);
        }
        events<<"{\"kind\":\"input.pointer\",\"frame\":"<<frame<<",\"command\":"<<next
              <<",\"op\":"<<quote(c.op)<<",\"stage\":"<<stage
              <<",\"mouse\":["<<injected_mouse.x<<','<<injected_mouse.y<<"]}\n";
        if(++stage>=6)complete_command=true;
        return;
    }
    // A real input sequence, observed hit rectangle -> down -> up. The
    // production Button/Selectable/canvas path remains the sole state writer.
    if(stage==1)io.AddMouseButtonEvent(0,true);
    if(stage==2) {
        io.AddMouseButtonEvent(0,false);
    }
    if(c.op=="text") {
        if(stage==4){io.AddKeyEvent(ImGuiMod_Ctrl,true);io.AddKeyEvent(ImGuiKey_A,true);}
        if(stage==5){io.AddKeyEvent(ImGuiKey_A,false);io.AddKeyEvent(ImGuiMod_Ctrl,false);}
        if(stage==6)io.AddInputCharactersUTF8(c.value.c_str());
        if(stage==8)io.AddKeyEvent(ImGuiKey_Enter,true);
        if(stage==9)io.AddKeyEvent(ImGuiKey_Enter,false);
        if(++stage>=12)complete_command=true;
    } else if(++stage>=5)complete_command=true;
}
void evaluated(std::size_t mo,const char* reason,float milliseconds) {
    if(!enabled)return;
    ++generation;volume_mo=mo;evaluation_reason=reason;kernel_ms=milliseconds;
}
void orbital_identity(const std::string& set,const std::string& dataset,
                      const std::string& spin,std::size_t source_index,
                      const std::string& association,
                      const std::string& coefficient_source,
                      bool direct_fchk_coefficients,bool density_verified) {
    if(!enabled)return;
    active_set=set;active_dataset=dataset;active_spin=spin;
    active_source_index=source_index;active_association=association;
    active_coefficient_source=coefficient_source;
    active_direct_fchk=direct_fchk_coefficients;
    active_density_verified=density_verified;
}
void ui_frame(std::size_t drawn,std::size_t requested) {
    drawn_ui_mo=drawn;requested_mo=requested;
}
void after_scene(const VolumeRenderer& renderer,const GridBox& box,std::size_t mo,std::size_t field_index) {
    if(!enabled)return;
    rendered_mo=mo;rendered_generation=generation;
    rendered_set=active_set;rendered_dataset=active_dataset;
    rendered_spin=active_spin;rendered_source_index=active_source_index;
    rendered_association=active_association;
    rendered_coefficient_source=active_coefficient_source;
    rendered_direct_fchk=active_direct_fchk;
    rendered_density_verified=active_density_verified;
    if(done() || !volume_command(commands[next]) || stage<4)return;
    const auto& c=commands[next];
    if(!c.value.empty() && std::stoull(c.value)!=mo) {finish("failed","selected MO does not match requested readback");return;}
    const int nx=renderer.nx(),ny=renderer.ny(),nz=renderer.nz();
    std::vector<float> volume(static_cast<std::size_t>(nx)*ny*nz);
    GLint binding=0,alignment=0;
    glGetIntegerv(0x806A,&binding);glGetIntegerv(GL_PACK_ALIGNMENT,&alignment);
    glBindTexture(0x806F,renderer.volume_texture());glPixelStorei(GL_PACK_ALIGNMENT,1);
    // CudaOrbitalEvaluator::evaluate returned after CUDA resource unmap. This
    // is the same texture just consumed by render_volume, before later UI updates.
    glGetTexImage(0x806F,0,GL_RED,GL_FLOAT,volume.data());
    glBindTexture(0x806F,binding);glPixelStorei(GL_PACK_ALIGNMENT,alignment);
    if(glGetError()!=GL_NO_ERROR)throw std::runtime_error("actual renderer volume readback failed");
    // Full-grid evidence is an opt-in extension of the existing v1 plan. It
    // stores the same texture just rendered, including its original float bits.
    const bool full=c.op=="volume_full";
    const std::string field_suffix=field_index?"-field"+std::to_string(field_index):"";
    const std::string binary_name=c.id+field_suffix+".volume.f32";
    if(full)write_volume_binary(output/binary_name,volume);
    std::mt19937 rng(20260905);std::vector<std::uint32_t> indices;
    for(int i=0;i<8192;++i)indices.push_back(rng()%static_cast<std::uint32_t>(volume.size()));
    std::sort(indices.begin(),indices.end());indices.erase(std::unique(indices.begin(),indices.end()),indices.end());
    std::ofstream out(output/(c.id+field_suffix+".volume.json"));out<<std::setprecision(17);
    out<<"{\"schema\":1,\"frame\":"<<frame<<",\"generation\":"<<generation<<",\"rendered_mo\":"<<mo
       <<",\"field_index\":"<<field_index
       <<",\"orbital_set\":"<<quote(rendered_set)<<",\"dataset\":"<<quote(rendered_dataset)
       <<",\"spin\":"<<quote(rendered_spin)<<",\"source_index\":";
    if(rendered_source_index==std::numeric_limits<std::size_t>::max())out<<"null";else out<<rendered_source_index;
    out<<",\"association\":"<<quote(rendered_association)
       <<",\"coefficient_source\":"<<quote(rendered_coefficient_source)
       <<",\"direct_fchk_coefficients\":"<<(rendered_direct_fchk?"true":"false")
       <<",\"density_verified\":"<<(rendered_density_verified?"true":"false")
       <<",\"texture_id\":"<<renderer.volume_texture()<<",\"nx\":"<<nx<<",\"ny\":"<<ny<<",\"nz\":"<<nz
       <<",\"grid_box_bohr\":["<<box.min_x<<','<<box.min_y<<','<<box.min_z<<','<<box.max_x<<','<<box.max_y<<','<<box.max_z
       <<"],\"layout\":\"x fastest; coordinates use CUDA float interpolation i/(n-1)\"";
    if(full)out<<",\"full_grid\":{\"file\":"<<quote(binary_name)
               <<",\"scalar_type\":\"IEEE754-float32-little-endian\",\"point_count\":"<<volume.size()
               <<",\"byte_count\":"<<volume.size()*4<<",\"source\":\"actual renderer texture readback\"}";
    out<<",\"samples\":[";
    bool first=true;for(auto idx:indices){if(!first)out<<',';first=false;out<<'['<<idx<<','<<volume[idx]<<']';}out<<"]}";
    out.close();if(!out)throw std::runtime_error("texture evidence metadata write failed");
    complete_command=true;
}
void scene_view(const ViewerLayout& layout, const OrbitCamera& camera) {
    if (!enabled) return;
    GLint viewport[4]{};
    glGetIntegerv(GL_VIEWPORT, viewport);
    const auto projection = scene_projection(viewport[2], viewport[3], camera.fov_degrees);
    const auto& scene = layout.scene;
    targets["scene.viewport"] = {ImVec2(scene.x, scene.y), ImVec2(scene.x + scene.width, scene.y + scene.height),
        nullptr, ImRect(ImVec2(scene.x, scene.y), ImVec2(scene.x + scene.width, scene.y + scene.height))};
    std::ostringstream out; out << std::setprecision(17);
    out << "{\"viewport_source\":\"GL_VIEWPORT read after scene rendering\",\"framebuffer_size\":["
        << layout.framebuffer_width << ',' << layout.framebuffer_height << "],\"window_size\":["
        << layout.window_width << ',' << layout.window_height << "],\"gl_viewport_xywh\":["
        << viewport[0] << ',' << viewport[1] << ',' << viewport[2] << ',' << viewport[3]
        << "],\"logical_scene_xywh\":[" << scene.x << ',' << scene.y << ',' << scene.width << ',' << scene.height
        << "],\"control_panel_xywh\":[" << layout.controls.x << ',' << layout.controls.y << ','
        << layout.controls.width << ',' << layout.controls.height << "],\"camera\":{\"yaw\":" << camera.yaw
        << ",\"pitch\":" << camera.pitch << ",\"distance\":" << camera.distance
        << ",\"fov_degrees\":" << camera.fov_degrees << ",\"fov_scope\":\"shorter viewport dimension\"}"
        << ",\"projection_inputs\":{\"aspect\":" << projection.aspect << ",\"tan_half_vertical_fov\":"
        << projection.tan_half_vertical_fov << ",\"source\":\"shared renderer projection function; not uniform readback\"}}";
    scene_view_json = out.str();
}
void hit(const std::string& id,ImVec2 lo,ImVec2 hi) {
    if(!enabled)return;auto* w=ImGui::GetCurrentWindow();
    targets[id]={lo,hi,w,w->ClipRect};
}
void item(const std::string& id) {
    if(!enabled || ImGui::GetCurrentWindow()->SkipItems)return;
    hit(id,ImGui::GetItemRectMin(),ImGui::GetItemRectMax());
}
void chrome_hit(const std::string& id, ImVec2 lo, ImVec2 hi) {
    if(!enabled)return;
    auto* window=ImGui::GetCurrentWindow();
    targets[id]={lo,hi,window,window->OuterRectClipped};
}
void anchor(const std::string& id) {
    if(!enabled)return;const auto p=ImGui::GetCursorScreenPos();hit(id,p,ImVec2(p.x+20,p.y+4));
}
void record(const std::string& kind,const std::string& json) {
    if(enabled && forensic && (kind=="nbo.integration" || kind=="aomo.selection" ||
       kind=="chemistry.route" || kind=="input.density_evidence" || kind=="input.pi_topology_evidence"))return;
    if(enabled && (kind=="nbo.integration" || kind=="aomo.selection")){
        const auto name=(kind=="nbo.integration"?"integration-":"selection-")+std::to_string(frame)+".json";
        std::ofstream file(output/name);file<<json;
        if(!file)throw std::runtime_error("Cannot preserve integrated orbital evidence");
        events<<"{\"frame\":"<<frame<<",\"kind\":"<<quote(kind)<<",\"file\":"<<quote(name)<<"}\n";
        events.flush();return;
    }
    if(enabled && (kind=="input.density_evidence" || kind=="input.pi_topology_evidence")) {
        events<<"{\"frame\":"<<frame<<",\"kind\":"<<quote(kind)<<",\"data\":"<<json<<"}\n";
        events.flush();
        return; // Full input evidence belongs to the load event, not per-frame traces.
    }
    if(enabled)trace.push_back("{\"kind\":"+quote(kind)+",\"data\":"+json+"}");
    if(enabled && (kind=="export.actual" || kind=="nbo.export" ||
                   kind=="browser.copy" || kind=="details.copy" || kind=="aomo.copy" ||
                   kind=="nbo.attach" || kind=="nbo.attach.error" || kind=="scene.pick" ||
                   kind=="aomo.selection.error" || kind=="input.package.error" ||
                   kind=="input.numerical_diagnostics" || kind=="input.frame_tracking" ||
                   (kind=="diagram.cache" && json.find("false")!=std::string::npos))) {
        if(kind=="diagram.cache")++diagram_generation;
        events<<"{\"frame\":"<<frame<<",\"kind\":"<<quote(kind)<<",\"data\":"<<json<<"}\n";events.flush();
        if(kind=="export.actual" || kind=="nbo.export") {
            events<<"{\"frame\":"<<frame<<",\"kind\":\"export.frame-trace\",\"data\":[";
            bool first=true;
            for(const auto& entry:trace){if(!first)events<<',';first=false;events<<entry;}
            events<<"]}\n";events.flush();
        }
    }
}
void field(const std::string& label,const std::string& value) {
    if(enabled)record("draw.text","{\"label\":"+quote(label)+",\"value\":"+quote(value)+"}");
}
std::filesystem::path export_base(const std::filesystem::path& original) {
    return enabled?output/export_name:original;
}
void end_frame(int width,int height,std::size_t applied,const ui::OrbitalUIState& ui,const Wavefunction* wf) {
    if(!enabled)return;
    const auto state=state_json(applied,ui,wf);
    if(!forensic)frames<<state<<'\n';
    const auto write_inspection=[&](const std::string& id) {
        std::ofstream out(output/(id+".view.json"));
        out<<"{\"schema\":1,\"state\":"<<state<<",\"targets\":[";bool first=true;
        for(const auto& [name,t]:targets){if(!first)out<<',';first=false;out<<"{\"id\":"<<quote(name)
            <<",\"rect\":["<<t.lo.x<<','<<t.lo.y<<','<<t.hi.x<<','<<t.hi.y
            <<"],\"clip_rect\":["<<t.clip.Min.x<<','<<t.clip.Min.y<<','<<t.clip.Max.x<<','<<t.clip.Max.y
            <<"],\"visible\":"<<(point_visible(t)?"true":"false")<<'}';}
        out<<"],\"draw_trace\":[";first=true;
        for(const auto& x:trace){if(!first)out<<',';first=false;out<<x;}
        out<<"],\"rendered\":"<<capture_rendered_frame_json(ImGui::GetDrawData(),ImGui::GetIO().Fonts)<<'}';
        out.close();if(!out)throw std::runtime_error("cannot save final display inspection");
    };
    if(!done() && detail_page_pending) {
        const auto it=targets.find("diagram.details.window");
        if(it==targets.end() || !it->second.window){finish("failed","details disappeared during capture");return;}
        auto* w=it->second.window;
        std::ostringstream name;name<<commands[next].id<<"-p"<<std::setfill('0')<<std::setw(3)<<detail_page;
        write_inspection(name.str());
        detail_page_pending=false;
        if(w->Scroll.y>=w->ScrollMax.y-1)complete_command=true;
        else if(++detail_page>=64){finish("failed","details page limit reached");return;}
        else {
            detail_scroll_before=w->Scroll.y;stage=0;
            // Less than one visible page gives overlap between captured pages.
            const float step=std::max(1.0f,std::floor(w->InnerRect.GetHeight()/(5.0f*w->CalcFontSize())*0.65f));
            // input_frame clears backend events, so store the pending wheel.
            cooldown=static_cast<int>(step);
        }
    }
    if(!done() && complete_command) {
        const auto c=commands[next];
        if(c.op=="expand" && targets.contains(c.id+".closed")) {
            finish("failed","header remained closed after real pointer input");return;
        }
        if(c.op=="inspect" || (forensic && c.op=="capture"))write_inspection(c.id);
        if(c.op=="capture" || c.op=="volume_full") {
            framebuffer(output/(c.id+".bmp"),width,height);
            std::ofstream out(output/(c.id+".ui.json"));
            out<<"{\"state\":"<<state<<",\"targets\":[";bool first=true;
            for(const auto& [id,t]:targets){if(!first)out<<',';first=false;out<<"{\"id\":"<<quote(id)<<",\"rect\":["<<t.lo.x<<','<<t.lo.y<<','<<t.hi.x<<','<<t.hi.y
                <<"],\"clip_rect\":["<<t.clip.Min.x<<','<<t.clip.Min.y<<','<<t.clip.Max.x<<','<<t.clip.Max.y
                <<"],\"visible\":"<<(point_visible(t)?"true":"false")<<'}';}
            out<<"],\"draw_trace\":[";first=true;for(const auto& x:trace){if(!first)out<<',';first=false;out<<x;}out<<"]}";
        }
        finish("executed");
    }
    if(done()) {
        std::ofstream summary(output/"session.json");
        summary<<"{\"schema\":1,\"seconds\":"<<number(elapsed())<<",\"frames\":"<<frame
               <<",\"commands\":"<<commands.size()<<",\"failed_commands\":"<<failures
               <<",\"capture_status\":"<<quote(failures?"failed":"completed")<<",\"scientific_verdict\":\"external checker required\"}";
        frames.flush();
    }
}
}
