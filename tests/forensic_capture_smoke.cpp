#include "cov/forensic_capture.hpp"
#include <imgui_internal.h>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition,const char* message) {if(!condition)throw std::runtime_error(message);}
std::string glyph_record(const std::string& capture,unsigned cp) {
    const auto position=capture.find("\"codepoint\":"+std::to_string(cp)+",");
    if(position==std::string::npos)return {};
    const auto begin=capture.rfind('{',position),end=capture.find('}',position);
    return capture.substr(begin,end-begin+1);
}
bool font_supports(const char* path,const ImWchar* ranges,std::initializer_list<unsigned> required) {
    if(!std::ifstream(path,std::ios::binary).good())return false;
    ImFontAtlas probe;auto* font=probe.AddFontFromFileTTF(path,13,nullptr,ranges);
    if(!font||!probe.Build())return false;
    for(auto cp:required)if(!font->FindGlyphNoFallback(static_cast<ImWchar>(cp)))return false;
    return true;
}
const char* scientific_font(const ImWchar* ranges) {
    // Existence does not imply coverage: Arial on Windows lacks subscript two.
    for(const char* path:{"C:/Windows/Fonts/seguisym.ttf","C:/Windows/Fonts/segoeui.ttf","C:/Windows/Fonts/cambria.ttc","C:/Windows/Fonts/arial.ttf","/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"})
        if(font_supports(path,ranges,{0x03c0,0x03c3,0x2082}))return path;
    return nullptr;
}
const char* cjk_font(const ImWchar* ranges) {
    for(const char* path:{"C:/Windows/Fonts/msyh.ttc","C:/Windows/Fonts/msyh.ttf","C:/Windows/Fonts/simhei.ttf","/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc","/usr/share/fonts/opentype/noto/NotoSansCJKsc-Regular.otf"})
        if(font_supports(path,ranges,{0x4e2d,0x6587}))return path;
    return nullptr;
}
}
int main(int argc,char** argv) {
    try {
        ImGui::CreateContext();auto& io=ImGui::GetIO();io.IniFilename=nullptr;io.DisplaySize={800,600};io.DeltaTime=1.f/60;
        auto* normal=io.Fonts->AddFontDefault();
        // Exercise the production merged-font arrangement, plus a separate
        // smaller font for scientific subscripts. Use actual font glyph shapes.
        static const ImWchar science_ranges[]={0x20,0xff,0x03a3,0x03c9,0x2080,0x2089,0};
        const auto* font_path=scientific_font(science_ranges);require(font_path,"scientific smoke fixture needs a font with actual pi, sigma and subscript two glyphs");
        ImFontConfig merged;merged.MergeMode=true;
        require(io.Fonts->AddFontFromFileTTF(font_path,13,&merged,science_ranges)!=nullptr,"merge font load failed");
        static const ImWchar cjk_ranges[]={0x4e2d,0x4e2d,0x6587,0x6587,0};
        const auto* cjk_path=cjk_font(cjk_ranges);require(cjk_path,"CJK smoke fixture needs Microsoft YaHei, SimHei or Noto CJK with Chinese glyphs");
        require(io.Fonts->AddFontFromFileTTF(cjk_path,13,&merged,cjk_ranges)!=nullptr,"CJK merge font load failed");
        auto* small=io.Fonts->AddFontFromFileTTF(font_path,9,nullptr,science_ranges);require(small,"small font load failed");
        require(io.Fonts->Build(),"font atlas build failed");io.Fonts->SetTexID(reinterpret_cast<ImTextureID>(std::uintptr_t(1)));
        require(normal->FindGlyphNoFallback(0x03c0)&&small->FindGlyphNoFallback(0x2082),"science fixture font lacks pi/subscript two");
        require(normal->FindGlyphNoFallback(0x4e2d)&&normal->FindGlyphNoFallback(0x6587),"CJK fixture glyphs unavailable");
        ImGui::NewFrame();ImGui::SetNextWindowPos({10,10});ImGui::SetNextWindowSize({500,400});
        ImGui::Begin("render capture",nullptr,ImGuiWindowFlags_NoSavedSettings);
        ImGui::TextUnformatted("widget label");auto* draw=ImGui::GetWindowDrawList();
        draw->AddText(normal,13,{40,90},IM_COL32_WHITE,"raw graph");
        draw->AddText(normal,13,{40,120},IM_COL32_WHITE,"\xcf\x83\xcf\x80");
        draw->AddText(small,9,{56,126},IM_COL32_WHITE,"\xe2\x82\x82");
        draw->AddText(normal,13,{40,300},IM_COL32_WHITE,"\xe4\xb8\xad\xe6\x96\x87");
        // The source character is intentionally unavailable. Capture must
        // expose the rendered fallback, never the original input character.
        require(!normal->FindGlyphNoFallback(0x2603),"fallback fixture unexpectedly contains snowman");
        draw->AddText(normal,13,{40,150},IM_COL32_WHITE,"\xe2\x98\x83");
        // CPU fine clipping should identify the glyph but mark it partial.
        const ImVec4 fine_clip{42,180,44,193};
        draw->AddText(normal,13,{40,180},IM_COL32_WHITE,"Q",nullptr,0,&fine_clip);
        // GPU clip keeps a complete quad outside its scissor rectangle.
        draw->PushClipRect({40,220},{50,240},true);draw->AddText(normal,13,{46,224},IM_COL32_WHITE,"W");draw->PopClipRect();
        draw->AddRectFilled({40,260},{80,280},IM_COL32_WHITE);
        draw->PushTextureID(reinterpret_cast<ImTextureID>(std::uintptr_t(2)));draw->AddImage(reinterpret_cast<ImTextureID>(std::uintptr_t(2)),{90,260},{110,280});draw->PopTextureID();
        // A later translucent window overlaps graph text; this must remain
        // uncertain, even though no opaque-pixel coverage was proven.
        ImGui::End();ImGui::SetNextWindowPos({36,86});ImGui::SetNextWindowSize({160,55});ImGui::SetNextWindowBgAlpha(0.25f);
        ImGui::Begin("later translucent",nullptr,ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoTitleBar);ImGui::TextUnformatted("cover");ImGui::End();
        auto* foreground=ImGui::GetForegroundDrawList();foreground->PushClipRect({0,0},{1000,800},false);
        foreground->AddText(normal,13,{900,450},IM_COL32_WHITE,"Z");
        foreground->AddText(normal,13,{797,480},IM_COL32_WHITE,"V");foreground->PopClipRect();ImGui::Render();
        cov::validation::RenderedFrameCaptureLimits full_limits;full_limits.include_glyph_records=true;
        const auto json=cov::validation::capture_rendered_frame_json(ImGui::GetDrawData(),io.Fonts,full_limits);
        require(json.find("\"text\":\"widget label\"")!=std::string::npos,"widget label was not reconstructed");
        require(json.find("\"text\":\"raw graph\"")!=std::string::npos,"raw AddText graph label was not reconstructed");
        require(json.find("\"codepoint\":960,")!=std::string::npos,"merged pi glyph was not decoded");
        require(json.find("\"codepoint\":8322,")!=std::string::npos,"subscript glyph was not decoded");
        require(json.find("\"text\":\"\xe4\xb8\xad\xe6\x96\x87\"")!=std::string::npos,"actual merged Chinese glyphs were not reconstructed");
        require(json.find("\"codepoint\":20013,")!=std::string::npos&&json.find("\"codepoint\":25991,")!=std::string::npos,"Chinese codepoint evidence absent");
        require(json.find("\"codepoint\":9731,")==std::string::npos,"capture invented unavailable snowman");
        require(json.find("\"codepoint\":"+std::to_string(normal->FallbackGlyph->Codepoint)+",")!=std::string::npos,"actual displayed fallback glyph absent");
        require(glyph_record(json,'Q').find("\"clipping\":\"software_partial\"")!=std::string::npos,"fine clipped glyph advertised as whole");
        require(glyph_record(json,'W').find("\"clipping\":\"command_partial\"")!=std::string::npos,"scissor clipped glyph advertised as whole");
        require(glyph_record(json,'Z').find("\"clipping\":\"outside_viewport\"")!=std::string::npos,"outside viewport glyph advertised as visible");
        require(glyph_record(json,'V').find("\"clipping\":\"viewport_partial\"")!=std::string::npos,"viewport boundary clipping not distinguished");
        require(json.find("possible_later_window_cover")!=std::string::npos,"later translucent window uncertainty absent");
        const auto sub=glyph_record(json,8322),pi=glyph_record(json,960);
        require(sub.find("\"font_size\":9,")!=std::string::npos&&pi.find("\"font_size\":13,")!=std::string::npos,"scientific font sizes lost");
        require(json.find("\"nonfont_commands\":0")==std::string::npos,"nonfont texture command not excluded");
        const auto compact=cov::validation::capture_rendered_frame_json(ImGui::GetDrawData(),io.Fonts);
        require(compact.find("\"glyph_records_included\":false,\"glyphs\":[]")!=std::string::npos,"default capture unexpectedly includes per-glyph records");
        require(compact.find("\"text\":\"raw graph\"")!=std::string::npos,"default capture lost text runs");
        cov::validation::RenderedFrameCaptureLimits tiny;tiny.max_glyphs=1;
        require(cov::validation::capture_rendered_frame_json(ImGui::GetDrawData(),io.Fonts,tiny).find("\"status\":\"truncated\"")!=std::string::npos,"glyph limit not reported");
        require(cov::validation::capture_rendered_frame_json(nullptr,io.Fonts).find("\"status\":\"unavailable\"")!=std::string::npos,"null frame not reported");
        if(argc>1&&std::string(argv[1])=="--json")std::cout<<json<<'\n';
        else std::cout<<"forensic_capture_smoke passed (widgets, raw graph, merged science font, sizes, fallback, clipping, viewport, translucent cover, limits)\n";
        ImGui::DestroyContext();return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';if(ImGui::GetCurrentContext())ImGui::DestroyContext();return 1;}
}
