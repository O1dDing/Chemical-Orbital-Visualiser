#include "cov/forensic_capture.hpp"
#include <imgui_internal.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <locale>
#include <map>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace cov::validation {
namespace {
struct Rect { float x0=0,y0=0,x1=0,y1=0; };
Rect intersect(Rect a,Rect b) { return {std::max(a.x0,b.x0),std::max(a.y0,b.y0),std::min(a.x1,b.x1),std::min(a.y1,b.y1)}; }
bool positive(Rect r) { return r.x1>r.x0 && r.y1>r.y0; }
bool near(float a,float b,float eps=0.01f) { return std::abs(a-b)<=eps; }
bool same(Rect a,Rect b) { return near(a.x0,b.x0)&&near(a.y0,b.y0)&&near(a.x1,b.x1)&&near(a.y1,b.y1); }
std::string_view bounded_name(const char* name) {
    if(!name)return {};
    std::string_view s(name);std::size_t n=std::min<std::size_t>(512,s.size());
    if(n<s.size())while(n>0&&(static_cast<unsigned char>(s[n])&0xc0)==0x80)--n;
    return s.substr(0,n);
}
Rect rect(const ImRect& r) { return {r.Min.x,r.Min.y,r.Max.x,r.Max.y}; }
void box(std::ostream& o,Rect r) { o<<'['<<r.x0<<','<<r.y0<<','<<r.x1<<','<<r.y1<<']'; }
void quoted(std::ostream& o,std::string_view s) {
    o<<'"';
    for (unsigned char c:s) {
        switch(c) { case '"':o<<"\\\"";break;case '\\':o<<"\\\\";break;
        case '\n':o<<"\\n";break;case '\r':o<<"\\r";break;case '\t':o<<"\\t";break;
        default:if(c<32) o<<"\\u00"<<"0123456789abcdef"[c>>4]<<"0123456789abcdef"[c&15];else o<<char(c); }
    }
    o<<'"';
}
std::string utf8(unsigned cp) {
    std::string s;
    if(cp<=0x7f)s+=char(cp);
    else if(cp<=0x7ff){s+=char(0xc0|(cp>>6));s+=char(0x80|(cp&63));}
    else if(cp<=0xffff){s+=char(0xe0|(cp>>12));s+=char(0x80|((cp>>6)&63));s+=char(0x80|(cp&63));}
    else if(cp<=0x10ffff){s+=char(0xf0|(cp>>18));s+=char(0x80|((cp>>12)&63));s+=char(0x80|((cp>>6)&63));s+=char(0x80|(cp&63));}
    return s;
}
struct AtlasGlyph { const ImFontGlyph* glyph;const ImFont* font;int font_index;Rect uv; };
using UVKey=std::array<long long,4>;
UVKey key(Rect uv,const ImFontAtlas& a) { return {std::llround(uv.x0*a.TexWidth*64),std::llround(uv.y0*a.TexHeight*64),std::llround(uv.x1*a.TexWidth*64),std::llround(uv.y1*a.TexHeight*64)}; }
struct Window { const ImGuiWindow* window;int stack,draw_list; };
struct Glyph {
    int list=0,command=0;unsigned index_offset=0,codepoint=0;
    Rect bounds,full_bounds,clip,uv;
    float size=0,baseline=0,origin=0,advance=0,space=0;
    unsigned alpha=0;int font=-1;
    bool software_partial=false,ambiguous=false;
    std::string text,clipping,overlap;std::vector<int> covers;
};
struct Run { std::string text,clipping,overlap;int list=0,command=0,font=-1;float size=0,baseline=0;Rect bounds,clip;unsigned index_offset=0;std::size_t first=0,count=0;bool inferred_spaces=false,ambiguous=false; };
std::ostringstream stream() {std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(7);return s;}
}

std::string capture_rendered_frame_json(const ImDrawData* data,ImFontAtlas* atlas) {
    return capture_rendered_frame_json(data,atlas,{});
}
std::string capture_rendered_frame_json(const ImDrawData* data,ImFontAtlas* atlas,const RenderedFrameCaptureLimits& limits) {
    if(!data||!data->Valid||!atlas||!atlas->IsBuilt()||atlas->TexWidth<=0||atlas->TexHeight<=0)
        return "{\"schema\":\"cov.rendered-frame.v1\",\"status\":\"unavailable\",\"reason\":\"valid_rendered_draw_data_and_built_atlas_required\",\"glyphs\":[],\"text_runs\":[],\"windows\":[]}";
    const Rect viewport{data->DisplayPos.x,data->DisplayPos.y,data->DisplayPos.x+data->DisplaySize.x,data->DisplayPos.y+data->DisplaySize.y};
    std::vector<AtlasGlyph> dictionary;std::map<UVKey,std::vector<std::size_t>> exact;
    for(int f=0;f<atlas->Fonts.Size;++f) {
        const auto* font=atlas->Fonts[f];
        for(const auto& g:font->Glyphs) if(g.Visible&&g.U1>g.U0&&g.V1>g.V0) {
            const Rect uv{g.U0,g.V0,g.U1,g.V1};exact[key(uv,*atlas)].push_back(dictionary.size());dictionary.push_back({&g,font,f,uv});
        }
    }
    std::unordered_map<const ImDrawList*,int> list_indices;
    for(int l=0;l<data->CmdListsCount;++l) list_indices[data->CmdLists[l]]=l;
    std::vector<Window> windows;bool truncated=false,windows_truncated=false;
    const auto* context=ImGui::GetCurrentContext();
    if(context) for(int w=0;w<context->Windows.Size;++w) {
        if(windows.size()>=limits.max_windows){truncated=true;windows_truncated=true;break;}
        auto* win=context->Windows[w];auto found=list_indices.find(win->DrawList);
        windows.push_back({win,w,found==list_indices.end()?-1:found->second});
    }
    std::vector<Glyph> glyphs;std::size_t scanned=0,nonfont=0,unknown=0,ambiguous=0,callbacks=0,invalid=0,text_bytes=0;
    bool stop=false;
    for(int l=0;l<data->CmdListsCount&&!stop;++l) {
        const auto* list=data->CmdLists[l];
        for(int c=0;c<list->CmdBuffer.Size&&!stop;++c) {
            const auto& cmd=list->CmdBuffer[c];
            if(cmd.UserCallback){++callbacks;continue;}
            if(cmd.GetTexID()!=atlas->TexID){++nonfont;continue;}
            const Rect clip{cmd.ClipRect.x,cmd.ClipRect.y,cmd.ClipRect.z,cmd.ClipRect.w};
            unsigned i=0;
            while(i+5<cmd.ElemCount) {
                if(scanned+6>limits.max_index_elements||glyphs.size()>=limits.max_glyphs){truncated=true;stop=true;break;}
                scanned+=6;
                const std::uint64_t begin=std::uint64_t(cmd.IdxOffset)+i;
                if(begin+5>=std::uint64_t(list->IdxBuffer.Size)){++invalid;break;}
                const auto* ids=list->IdxBuffer.Data+begin;
                // ImGui AddText emits two triangles sharing vertices 0 and 2.
                if(ids[0]!=ids[3]||ids[2]!=ids[4]||ids[0]==ids[1]||ids[0]==ids[2]||ids[0]==ids[5]||ids[1]==ids[2]||ids[1]==ids[5]||ids[2]==ids[5]) {i+=3;continue;}
                std::array<const ImDrawVert*,4> v{};const std::array<unsigned,4> offset{0,1,2,5};bool valid=true;
                for(int k=0;k<4;++k) {
                    const std::uint64_t index=std::uint64_t(cmd.VtxOffset)+ids[offset[k]];
                    if(index>=std::uint64_t(list->VtxBuffer.Size)){valid=false;break;}v[k]=&list->VtxBuffer[int(index)];
                }
                if(!valid){++invalid;i+=6;continue;}
                const Rect b{v[0]->pos.x,v[0]->pos.y,v[2]->pos.x,v[2]->pos.y},uv{v[0]->uv.x,v[0]->uv.y,v[2]->uv.x,v[2]->uv.y};
                bool quad=positive(b)&&positive(uv)&&near(v[1]->pos.x,b.x1)&&near(v[1]->pos.y,b.y0)&&near(v[3]->pos.x,b.x0)&&near(v[3]->pos.y,b.y1)&&near(v[1]->uv.x,uv.x1,1e-6f)&&near(v[1]->uv.y,uv.y0,1e-6f)&&near(v[3]->uv.x,uv.x0,1e-6f)&&near(v[3]->uv.y,uv.y1,1e-6f);
                for(int k=1;k<4;++k)quad=quad&&v[k]->col==v[0]->col;
                if(!quad){i+=3;continue;}
                i+=6;const unsigned alpha=(v[0]->col>>IM_COL32_A_SHIFT)&255;
                if(!alpha)continue;
                std::vector<std::size_t> candidates;auto found=exact.find(key(uv,*atlas));
                if(found!=exact.end())candidates=found->second;
                bool partial=false;
                if(candidates.empty()) {
                    // CPU fine clipping changes UVs as well as positions. Only
                    // infer a glyph if its atlas rectangle contains this UV box.
                    const float eps_u=0.01f/atlas->TexWidth,eps_v=0.01f/atlas->TexHeight;
                    for(std::size_t d=0;d<dictionary.size();++d) {
                        const auto& g=dictionary[d];const auto* gg=g.glyph;
                        if(uv.x0<g.uv.x0-eps_u||uv.y0<g.uv.y0-eps_v||uv.x1>g.uv.x1+eps_u||uv.y1>g.uv.y1+eps_v)continue;
                        const float sx=(b.x1-b.x0)/(uv.x1-uv.x0)*(g.uv.x1-g.uv.x0)/(gg->X1-gg->X0);
                        const float sy=(b.y1-b.y0)/(uv.y1-uv.y0)*(g.uv.y1-g.uv.y0)/(gg->Y1-gg->Y0);
                        if(std::isfinite(sx)&&sx>0&&near(sx,sy,std::max(0.01f,sx*0.01f)))candidates.push_back(d);
                    }
                    partial=!candidates.empty();
                }
                if(candidates.empty()){++unknown;continue;}
                const auto& match=dictionary[candidates.front()];const auto& g=*match.glyph;
                const float scale=(b.x1-b.x0)/(uv.x1-uv.x0)*(match.uv.x1-match.uv.x0)/(g.X1-g.X0);
                const float scale_y=(b.y1-b.y0)/(uv.y1-uv.y0)*(match.uv.y1-match.uv.y0)/(g.Y1-g.Y0);
                if(!std::isfinite(scale)||scale<=0||!near(scale,scale_y,std::max(0.01f,scale*0.01f))){++unknown;continue;}
                const float full_x=b.x0-(uv.x0-match.uv.x0)/(match.uv.x1-match.uv.x0)*(g.X1-g.X0)*scale;
                const float full_y=b.y0-(uv.y0-match.uv.y0)/(match.uv.y1-match.uv.y0)*(g.Y1-g.Y0)*scale;
                Glyph out;out.list=l;out.command=c;out.index_offset=cmd.IdxOffset+i-6;out.codepoint=g.Codepoint;out.bounds=b;out.full_bounds={full_x,full_y,full_x+(g.X1-g.X0)*scale,full_y+(g.Y1-g.Y0)*scale};out.clip=clip;out.uv=uv;out.size=match.font->FontSize*scale;out.font=match.font_index;out.alpha=alpha;out.baseline=full_y-g.Y0*scale;out.origin=full_x-g.X0*scale;out.advance=g.AdvanceX*scale;out.space=match.font->GetCharAdvance(' ')*scale;out.software_partial=partial;
                for(auto d:candidates) {
                    const auto& other=*dictionary[d].glyph;
                    if(other.Codepoint!=g.Codepoint||!near(dictionary[d].font->FontSize,match.font->FontSize)||!near(other.X0,g.X0)||!near(other.Y0,g.Y0)||!near(other.AdvanceX,g.AdvanceX))out.ambiguous=true;
                }
                if(out.ambiguous){++ambiguous;out.text="\xef\xbf\xbd";}else out.text=utf8(g.Codepoint);
                if(text_bytes+out.text.size()>limits.max_text_bytes){truncated=true;stop=true;break;}text_bytes+=out.text.size();
                const auto command_visible=intersect(b,clip),visible=intersect(command_visible,viewport);
                if(!positive(visible))out.clipping=positive(command_visible)?"outside_viewport":"outside_command_clip";
                else if(partial)out.clipping="software_partial";
                else if(!same(visible,b))out.clipping=!same(command_visible,b)?"command_partial":"viewport_partial";
                else out.clipping="unclipped";
                out.overlap=context?(windows_truncated?"window_stack_truncated":"no_later_window_overlap"):"window_stack_unavailable";
                for(const auto& w:windows)if(w.draw_list>l&&w.window->Active&&!w.window->Hidden&&positive(intersect(visible,rect(w.window->OuterRectClipped))))out.covers.push_back(w.stack);
                if(!out.covers.empty())out.overlap="possible_later_window_cover";
                glyphs.push_back(std::move(out));
            }
        }
    }
    std::vector<Run> runs;
    for(std::size_t g=0;g<glyphs.size();++g) {
        const auto& glyph=glyphs[g];bool append=false;
        if(!runs.empty()&&g>0) {
            const auto& prev=glyphs[g-1];const auto& run=runs.back();const float gap=glyph.origin-(prev.origin+prev.advance);
            append=run.list==glyph.list&&run.command==glyph.command&&run.clipping==glyph.clipping&&run.overlap==glyph.overlap&&prev.font==glyph.font&&!prev.ambiguous&&!glyph.ambiguous&&near(run.size,glyph.size,0.05f)&&near(run.baseline,glyph.baseline,0.75f)&&gap>=-1.0f&&gap<std::max(glyph.size*3,prev.space*8);
            if(append&&gap>prev.space*0.55f&&prev.space>0) {const int spaces=std::clamp(int(std::lround(gap/prev.space)),1,8);runs.back().text.append(spaces,' ');runs.back().inferred_spaces=true;}
        }
        if(!append) {
            if(runs.size()>=limits.max_runs){truncated=true;break;}
            Run run;run.list=glyph.list;run.command=glyph.command;run.font=glyph.font;run.size=glyph.size;run.baseline=glyph.baseline;run.bounds=glyph.bounds;run.clip=glyph.clip;run.index_offset=glyph.index_offset;run.clipping=glyph.clipping;run.overlap=glyph.overlap;run.first=g;run.ambiguous=glyph.ambiguous;runs.push_back(std::move(run));
        }
        auto& run=runs.back();run.text+=glyph.text;++run.count;run.bounds={std::min(run.bounds.x0,glyph.bounds.x0),std::min(run.bounds.y0,glyph.bounds.y0),std::max(run.bounds.x1,glyph.bounds.x1),std::max(run.bounds.y1,glyph.bounds.y1)};
    }
    auto o=stream();o<<"{\"schema\":\"cov.rendered-frame.v1\",\"status\":";quoted(o,truncated?"truncated":"captured");
    o<<",\"source\":\"final_ImDrawData_font_quads\",\"legibility\":\"raster_unverified\",\"limitations\":[\"spaces_and_runs_inferred\",\"unrendered_characters_unrecoverable\",\"later_window_overlap_is_conservative_including_translucent_windows\",\"same_list_primitives_external_scene_and_callbacks_not_raster_resolved\"],\"viewport\":";box(o,viewport);
    o<<",\"framebuffer_scale\":["<<data->FramebufferScale.x<<','<<data->FramebufferScale.y<<"],\"draw_list_count\":"<<data->CmdListsCount<<",\"window_stack_available\":"<<(context?"true":"false")<<",\"counters\":{\"captured_glyphs\":"<<glyphs.size()<<",\"text_runs\":"<<runs.size()<<",\"scanned_index_elements\":"<<scanned<<",\"nonfont_commands\":"<<nonfont<<",\"unmatched_font_texture_quads\":"<<unknown<<",\"ambiguous_glyphs\":"<<ambiguous<<",\"callbacks_not_decoded\":"<<callbacks<<",\"invalid_index_groups\":"<<invalid<<"},\"windows\":[";
    bool comma=false;for(const auto& w:windows) {
        if(comma)o<<',';comma=true;const auto* win=w.window;o<<"{\"stack_order\":"<<w.stack<<",\"name\":";quoted(o,bounded_name(win->Name));o<<",\"draw_list\":"<<w.draw_list<<",\"active\":"<<(win->Active?"true":"false")<<",\"hidden\":"<<(win->Hidden?"true":"false")<<",\"collapsed\":"<<(win->Collapsed?"true":"false")<<",\"no_background\":"<<((win->Flags&ImGuiWindowFlags_NoBackground)?"true":"false")<<",\"rect\":";box(o,rect(win->Rect()));o<<",\"outer_clipped_rect\":";box(o,rect(win->OuterRectClipped));o<<",\"content_clip\":";box(o,rect(win->ClipRect));o<<",\"parent\":";quoted(o,bounded_name(win->ParentWindow?win->ParentWindow->Name:nullptr));o<<'}';
    }
    o<<"],\"glyph_records_included\":"<<(limits.include_glyph_records?"true":"false")<<",\"glyphs\":[";comma=false;
    if(limits.include_glyph_records)for(const auto& g:glyphs) {
        if(comma)o<<',';comma=true;o<<"{\"text\":";quoted(o,g.text);o<<",\"codepoint\":";if(g.ambiguous)o<<"null";else o<<g.codepoint;o<<",\"font_index\":"<<g.font<<",\"font_size\":"<<g.size<<",\"baseline_origin_y\":"<<g.baseline<<",\"draw_list\":"<<g.list<<",\"draw_command\":"<<g.command<<",\"index_offset\":"<<g.index_offset<<",\"owner\":";quoted(o,bounded_name(data->CmdLists[g.list]->_OwnerName));o<<",\"bounds\":";box(o,g.bounds);o<<",\"unclipped_glyph_bounds\":";box(o,g.full_bounds);o<<",\"clip\":";box(o,g.clip);o<<",\"uv\":";box(o,g.uv);o<<",\"alpha\":"<<g.alpha<<",\"ambiguous\":"<<(g.ambiguous?"true":"false")<<",\"clipping\":";quoted(o,g.clipping);o<<",\"occlusion\":";quoted(o,g.overlap);o<<",\"later_window_stack_orders\":[";for(std::size_t k=0;k<g.covers.size();++k){if(k)o<<',';o<<g.covers[k];}o<<"]}";
    }
    o<<"],\"text_runs\":[";comma=false;
    for(const auto& r:runs) {
        if(comma)o<<',';comma=true;o<<"{\"text\":";quoted(o,r.text);o<<",\"draw_list\":"<<r.list<<",\"draw_command\":"<<r.command<<",\"font_index\":"<<r.font<<",\"ambiguous\":"<<(r.ambiguous?"true":"false")<<",\"font_size\":"<<r.size<<",\"baseline_origin_y\":"<<r.baseline<<",\"bounds\":";box(o,r.bounds);o<<",\"clip\":";box(o,r.clip);o<<",\"index_offset\":"<<r.index_offset<<",\"owner\":";quoted(o,bounded_name(data->CmdLists[r.list]->_OwnerName));o<<",\"first_glyph\":"<<r.first<<",\"glyph_count\":"<<r.count<<",\"spaces_inferred\":"<<(r.inferred_spaces?"true":"false")<<",\"clipping\":";quoted(o,r.clipping);o<<",\"occlusion\":";quoted(o,r.overlap);o<<'}';
    }
    o<<"]}";auto result=o.str();
    if(result.size()>limits.max_json_bytes)return "{\"schema\":\"cov.rendered-frame.v1\",\"status\":\"truncated\",\"reason\":\"json_byte_limit_exceeded_reduce_capture_limits\",\"glyphs\":[],\"text_runs\":[],\"windows\":[]}";
    return result;
}
} // namespace cov::validation
