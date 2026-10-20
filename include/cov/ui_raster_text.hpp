#pragma once

#include <imgui.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include "cov/ui_scientific_text.hpp"

namespace cov::ui {

inline int raster_text_width(const std::string& text, int scale) {
    if(!ImGui::GetCurrentContext() || !ImGui::GetFont()) return 0;
    return static_cast<int>(std::ceil(scientific_text_layout(text,12.0f*scale).width));
}

// Use the viewer's loaded scientific and CJK glyphs in exported raster images.
inline void raster_plain_text(std::vector<std::uint8_t>& rgba,int width,int height,
                        float x,float y,const std::string& text,
                        std::uint8_t r,std::uint8_t g,std::uint8_t b,float size,
                        std::size_t* missing=nullptr) {
    if(!ImGui::GetCurrentContext() || !ImGui::GetFont()) return;
    unsigned char* atlas=nullptr; int aw=0,ah=0;
    ImGui::GetIO().Fonts->GetTexDataAsAlpha8(&atlas,&aw,&ah);
    const ImFont* font=ImGui::GetFont();
    if(!atlas || aw<=0 || ah<=0) return;
    const float scale=size/std::max(1.0f,font->FontSize);
    float cursor=x;
    for(std::size_t p=0;p<text.size();) {
        const auto c=static_cast<unsigned char>(text[p]);
        unsigned cp=c; std::size_t length=1;
        if((c&0xe0)==0xc0){cp=c&0x1f;length=2;}
        else if((c&0xf0)==0xe0){cp=c&0x0f;length=3;}
        else if((c&0xf8)==0xf0){cp=c&0x07;length=4;}
        if(p+length>text.size()){cp='?';length=1;}
        else for(std::size_t k=1;k<length;++k)
            cp=(cp<<6)|(static_cast<unsigned char>(text[p+k])&0x3f);
        p+=length;
        const auto* glyph=font->FindGlyphNoFallback(static_cast<ImWchar>(cp));
        if(!glyph){if(missing)++*missing;glyph=font->FindGlyph('?');}
        if(!glyph)continue;
        const int x0=static_cast<int>(std::floor(cursor+glyph->X0*scale));
        const int x1=static_cast<int>(std::ceil(cursor+glyph->X1*scale));
        const int y0=static_cast<int>(std::floor(y+glyph->Y0*scale));
        const int y1=static_cast<int>(std::ceil(y+glyph->Y1*scale));
        const int u0=static_cast<int>(glyph->U0*aw),u1=static_cast<int>(glyph->U1*aw);
        const int v0=static_cast<int>(glyph->V0*ah),v1=static_cast<int>(glyph->V1*ah);
        for(int yy=std::max(0,y0);yy<std::min(height,y1);++yy)
            for(int xx=std::max(0,x0);xx<std::min(width,x1);++xx) {
                const int u=std::clamp(u0+(xx-x0)*(u1-u0)/std::max(1,x1-x0),0,aw-1);
                const int v=std::clamp(v0+(yy-y0)*(v1-v0)/std::max(1,y1-y0),0,ah-1);
                const float a=atlas[v*aw+u]/255.0f;
                const auto at=(static_cast<std::size_t>(yy)*width+xx)*4;
                rgba[at]=static_cast<std::uint8_t>(rgba[at]*(1-a)+r*a);
                rgba[at+1]=static_cast<std::uint8_t>(rgba[at+1]*(1-a)+g*a);
                rgba[at+2]=static_cast<std::uint8_t>(rgba[at+2]*(1-a)+b*a);
            }
        cursor+=glyph->AdvanceX*scale;
    }
}
inline void raster_scientific_text(std::vector<std::uint8_t>& rgba,int width,int height,
                        float x,float y,const std::string& text,
                        std::uint8_t r,std::uint8_t g,std::uint8_t b,float size,
                        std::size_t* missing=nullptr,bool atom_indices=false) {
    if(!ImGui::GetCurrentContext() || !ImGui::GetFont())return;
    for(const auto& run:scientific_text_layout(text,size,false,atom_indices).runs)
        raster_plain_text(rgba,width,height,x+run.x,y+run.y,run.text,
                          r,g,b,size*run.scale,missing);
}
inline void raster_text(std::vector<std::uint8_t>& rgba,int width,int height,
                        int x,int y,const std::string& text,
                        std::uint8_t r,std::uint8_t g,std::uint8_t b,int size) {
    raster_scientific_text(rgba,width,height,static_cast<float>(x),
        static_cast<float>(y),text,r,g,b,12.0f*size);
}
} // namespace cov::ui
