/*
 * SPDX-License-Identifier: CPAL-1.0
 *
 * The contents of this file are subject to the Common Public Attribution License
 * Version 1.0 (the “License”); you may not use this file except in compliance with the
 * License. You may obtain a copy of the License at
 * https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
 * License Version 1.1 but Sections 14 and 15 have been added to cover use of software
 * over a computer network and provide for limited attribution for the Original
 * Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
 *
 * Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
 * WARRANTY OF ANY KIND, either express or implied. See the License for the specific
 * language governing rights and limitations under the License.
 *
 * The Original Code is Chemical Orbital Visualiser (COV).
 *
 * The Original Developer is O1dDing, who is also the Initial Developer.
 *
 * The Initial Developer of the Original Code is O1dDing. All portions of the code
 * written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
 *
 * Contributor O1dDing.
 */

#include "cov/scientific_text.hpp"
#include "cov/ui_scientific_text.hpp"
#include "cov/ui.hpp"
#include "cov/nbo_aomo_ui.hpp"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool ok,const char* why) { if(!ok)throw std::runtime_error(why); }
bool near(float a,float b) { return std::abs(a-b)<.02f; }
void token(const std::string& text,const std::string& base,
           const std::string& sub={},const std::string& sup={}) {
    const auto parts=cov::scientific_token(text);
    if(parts.size()!=1 || parts[0].base!=base || parts[0].subscript!=sub ||
       parts[0].superscript!=sup)throw std::runtime_error("Unexpected notation for "+text);
}
void plain(const std::string& text) {
    std::string restored;
    for(const auto& part:cov::scientific_text_parts(text)) {
        require(part.subscript.empty()&&part.superscript.empty(),"Prose or source identity gained a script");
        restored+=part.base;
    }
    require(restored==text,"Prose or source identity changed");
}
std::string read_file(const std::filesystem::path& path) {
    std::ifstream in(path,std::ios::binary);
    require(static_cast<bool>(in),"Unable to read exported typography artifact");
    return {std::istreambuf_iterator<char>(in),{}};
}
std::string svg_node(const std::string& svg,const std::string& id) {
    const auto start=svg.find("<g id=\""+id+"\"");
    require(start!=std::string::npos,"Export omitted a typography node");
    const auto end=svg.find("</g>",start);
    require(end!=std::string::npos,"Malformed exported typography node");
    return svg.substr(start,end-start);
}
void export_typography(const std::filesystem::path& base) {
    using namespace cov::ui;
    NboAomoViewSnapshot view;
    cov::NboIntegration integration;integration.id="synthetic-typography-only";
    view.id="synthetic-typography-only";view.integration_id=integration.id;
    view.capability_status="available";view.canvas_width=1100;view.canvas_height=620;
    view.footer_y=600;view.label_font_size=28;view.paper_export=true;
    view.captions.push_back({"title","Scientific notation: production SVG and PNG export",24,15,1000,30});
    view.captions.push_back({"atoms","Atomic identity and AO direction",35,65,480,30});
    view.captions.push_back({"irreps","Irrep, approximation and copy count",555,65,500,30});
    const auto add=[&](const std::string& id,const std::string& label,bool atoms,int column,int row) {
        NboAomoNode node;node.id=id;node.label=node.individual_label=label;
        node.label_atom_indices=atoms;node.lane=NboAomoLane::Centre;
        node.quantitative_energy=true;node.x=35+column*520.f;node.y=112+row*78.f;
        node.width=100;node.height=10;node.label_x=node.x;node.label_y=node.y+18;
        node.label_width=scientific_text_layout(label,view.label_font_size,false,atoms).width;
        node.label_height=view.label_font_size*1.08f;view.nodes.push_back(std::move(node));
    };
    add("typo-atom-H1","H1 2px / NAO 1",true,0,0);
    add("typo-atom-B1","B1 2px / NAO 2",true,0,1);
    add("typo-atom-F1","F1 2px / NAO 3",true,0,2);
    add("typo-atom-d","Cr1 d_x2-y2",true,0,3);
    add("typo-atom-poly","Gaussian AO 8 / shell 3 d x^2",true,0,4);
    add("typo-irrep-B1","B1",false,1,0);
    add("typo-irrep-numbered","12t₂g [alpha]",false,1,1);
    add("typo-irrep-approx","≈t₂g [MO 8]",false,1,2);
    add("typo-irrep-linear","≈σg+ [MO 9]",false,1,3);
    add("typo-irrep-copy","2×t₂g (6)",false,1,4);
    add("typo-irrep-prime","E2''",false,1,5);
    const auto result=export_nbo_aomo_bundle(view,integration,base,cov::DiagramExportContent::Images);
    if(!result.svg||!result.png)throw std::runtime_error("Production typography export failed: "+result.error);
    require(!result.json&&!result.csv,"Images-only typography fixture exported analysis data");
    const auto svg=read_file(result.svg_path);
    for(const std::string atom:{"H1","B1","F1"}) {
        const auto node=svg_node(svg,"typo-atom-"+atom);
        require(node.find(">"+atom+"</tspan>")!=std::string::npos&&
            node.find(">2p</tspan>")!=std::string::npos&&node.find(">x</tspan>")!=std::string::npos,
            "Production SVG lost atomic baseline context or AO direction subscript");
    }
    const auto irrep=svg_node(svg,"typo-irrep-B1");
    require(irrep.find(">B</tspan>")!=std::string::npos&&
        irrep.find(">1</tspan>")!=std::string::npos&&irrep.find(">B1</tspan>")==std::string::npos,
        "Production SVG confused a pure irrep with an atom identity");
    for(const std::string id:{"typo-irrep-approx","typo-irrep-copy"})
        require(svg_node(svg,id).find(">2g</tspan>")!=std::string::npos,
            "Production SVG lost the full subscript after a baseline prefix");
    const auto linear=svg_node(svg,"typo-irrep-linear");
    require(linear.find(">g</tspan>")!=std::string::npos&&linear.find(">+</tspan>")!=std::string::npos,
        "Production SVG lost prefixed linear parity or reflection script");
    const auto png=read_file(result.png_path);
    require(png.size()>24&&png.compare(0,8,std::string("\x89PNG\r\n\x1a\n",8))==0,
        "Production typography raster is not a PNG");
    std::cout<<"production_preview="<<result.png_path.string()<<"; svg="<<result.svg_path.string()<<'\n';
}
struct Bounds { float left=1e6f,top=1e6f,right=-1e6f,bottom=-1e6f; };
Bounds bounds(const ImDrawList& draw,int first,int last) {
    Bounds result;
    for(int i=first;i<last;++i) {
        const auto p=draw.VtxBuffer[i].pos;
        result.left=std::min(result.left,p.x);result.top=std::min(result.top,p.y);
        result.right=std::max(result.right,p.x);result.bottom=std::max(result.bottom,p.y);
    }
    return result;
}

// Text creates axis-aligned, atlas-textured quads. Rasterizing those actual
// ImGui vertices checks the same positioned glyph geometry used by the viewer.
void preview(const ImDrawList& draw,const std::filesystem::path& path,int width,int height) {
    unsigned char* atlas=nullptr;int aw=0,ah=0;
    ImGui::GetIO().Fonts->GetTexDataAsAlpha8(&atlas,&aw,&ah);
    require(atlas&&aw>0&&ah>0,"Missing font atlas");
    std::vector<unsigned char> pixels(static_cast<std::size_t>(width)*height*3,255);
    for(int i=0;i+3<draw.VtxBuffer.Size;i+=4) {
        const auto& a=draw.VtxBuffer[i];const auto& c=draw.VtxBuffer[i+2];
        if(c.pos.x<=a.pos.x||c.pos.y<=a.pos.y)continue;
        const unsigned red=(a.col>>IM_COL32_R_SHIFT)&255;
        const unsigned green=(a.col>>IM_COL32_G_SHIFT)&255;
        const unsigned blue=(a.col>>IM_COL32_B_SHIFT)&255;
        const float opacity=((a.col>>IM_COL32_A_SHIFT)&255)/255.f;
        for(int y=std::max(0,static_cast<int>(std::floor(a.pos.y)));
            y<std::min(height,static_cast<int>(std::ceil(c.pos.y)));++y)
            for(int x=std::max(0,static_cast<int>(std::floor(a.pos.x)));
                x<std::min(width,static_cast<int>(std::ceil(c.pos.x)));++x) {
                const float u=a.uv.x+(c.uv.x-a.uv.x)*(x+.5f-a.pos.x)/(c.pos.x-a.pos.x);
                const float v=a.uv.y+(c.uv.y-a.uv.y)*(y+.5f-a.pos.y)/(c.pos.y-a.pos.y);
                const int tx=std::clamp(static_cast<int>(u*aw),0,aw-1);
                const int ty=std::clamp(static_cast<int>(v*ah),0,ah-1);
                const float alpha=opacity*atlas[ty*aw+tx]/255.f;
                const auto at=(static_cast<std::size_t>(y)*width+x)*3;
                pixels[at]=static_cast<unsigned char>(pixels[at]*(1-alpha)+red*alpha);
                pixels[at+1]=static_cast<unsigned char>(pixels[at+1]*(1-alpha)+green*alpha);
                pixels[at+2]=static_cast<unsigned char>(pixels[at+2]*(1-alpha)+blue*alpha);
            }
    }
    std::ofstream out(path,std::ios::binary);
    out<<"P6\n"<<width<<' '<<height<<"\n255\n";
    out.write(reinterpret_cast<const char*>(pixels.data()),static_cast<std::streamsize>(pixels.size()));
    require(static_cast<bool>(out),"Unable to write typography preview");
}
}

int main(int argc,char** argv) { try {
    token("12t₂g","12t","2g");token("12T2g","12T","2g");
    token("t₁u","t","1u");
    token("B1","B","1"); // The same spelling is an irrep outside atom context.
    for(const std::string atom:{"H1","B1","F1"}) {
        const auto parts=cov::scientific_text_parts(atom+" 2px / NAO 1",false,true);
        require(parts.size()>3&&parts[0].base==atom&&parts[0].subscript.empty()&&
            parts[0].superscript.empty()&&parts[2].base=="2p"&&parts[2].subscript=="x",
            "Typed atom context changed an atom index or suppressed the AO subscript");
    }
    const auto explicit_irrep=cov::scientific_token("B₁",false,true);
    require(explicit_irrep.size()==1&&explicit_irrep[0].base=="B"&&
        explicit_irrep[0].subscript=="1","Explicitly scripted irrep was suppressed in atom context");
    for(const std::string prefix:{"≈","2×"}) {
        const auto parts=cov::scientific_token(prefix+"t₂g");
        require(parts.size()==2&&parts[0].base==prefix&&parts[0].subscript.empty()&&
            parts[0].superscript.empty()&&parts[1].base=="t"&&parts[1].subscript=="2g",
            "Approximation or copy count prefix swallowed an irrep subscript");
    }
    const auto approximate_linear=cov::scientific_token("≈σg+");
    require(approximate_linear.size()==2&&approximate_linear[0].base=="≈"&&
        approximate_linear[1].base=="σ"&&approximate_linear[1].subscript=="g"&&
        approximate_linear[1].superscript=="+","Approximation prefix swallowed linear scripts");
    token("2σg+","2σ","g","+");token("Σ_u^-","Σ","u","−");
    token("πᵤ","π","u");token("σg⁺","σ","g","+");
    token("E2''","E","2","′′");token("a₁′","a","1","′");
    token("E″","E",{},"′′");token("t2g*","t","2g","*");
    token("C2v","C","2v");token("D3h","D","3h");
    token("Oh","O","h");token("D∞h","D","∞h");
    token("Cinfv","C","∞v");token("Cs","C","s");token("Ci","C","i");
    token("2px","2p","x");token("10p_x","10p","x");
    token("dxy","d","xy");token("d_x2-y2","d","x²−y²");
    token("dz2","d","z²");token("pₓ","p","x");
    const auto polynomial=cov::scientific_token("x^2y^3z");
    require(polynomial.size()==3&&polynomial[0].base=="x"&&polynomial[0].superscript=="2"&&
        polynomial[1].base=="y"&&polynomial[1].superscript=="3"&&
        polynomial[2].base=="z"&&polynomial[2].superscript.empty(),"Cartesian powers were not separated");
    plain("MO 8 [alpha] / NAO 12 / shell 3");
    plain("C1 C2 S4 / atom 12");
    const auto bare_group=cov::scientific_token("D2",true);
    require(bare_group.size()==1&&bare_group[0].base=="D"&&bare_group[0].subscript=="2",
        "Explicit point group was mistaken for an atom index");
    plain("source_id=GaussianAO:12; internal_index=7; atom=Cr1");
    plain("SO(3), d real m=1 cos; mixed directions; unknown prose");
    plain("A_g2"); // Unknown suffix order stays literal rather than being reordered.
    const std::string source="Cr1 dxy / NAO 1";
    const auto source_parts=cov::scientific_text_parts(source);
    require(source=="Cr1 dxy / NAO 1"&&source_parts[2].base=="d"&&
        source_parts[2].subscript=="xy","AO direction or source identity was lost");

    IMGUI_CHECKVERSION();ImGui::CreateContext();
    require(cov::ui::configure_fonts(24.f),"Unable to configure application fonts");
    auto& io=ImGui::GetIO();io.IniFilename=nullptr;io.DisplaySize={1100,400};io.DeltaTime=1.f/60;
    unsigned char* atlas=nullptr;int aw=0,ah=0;
    io.Fonts->GetTexDataAsAlpha8(&atlas,&aw,&ah);
    ImGui::NewFrame();
    const ImFont* font=ImGui::GetFont();
    for(ImWchar glyph:{0x03c3,0x03a3,0x221e,0x2032,0x2212,0x00b2,0x00b3})
        require(font->FindGlyphNoFallback(glyph)!=nullptr,"Required scientific glyph is unavailable");
    const auto finite=cov::ui::scientific_text_layout("12t₂g [alpha]",24);
    require(finite.runs.size()>=3&&finite.runs[0].text=="12t"&&
        finite.runs[0].y==0&&finite.runs[0].scale==1&&finite.runs[1].text=="2g"&&
        finite.runs[1].y>0&&finite.runs[1].scale<1,"Ordinal or full parity suffix has wrong placement");
    for(const std::string atom:{"H1","B1","F1"}) {
        const auto layout=cov::ui::scientific_text_layout(atom+" 2px",24,false,true);
        require(layout.runs.size()==4&&layout.runs[0].text==atom&&layout.runs[0].y==0&&
            layout.runs[0].scale==1&&layout.runs[2].text=="2p"&&layout.runs[2].y==0&&
            layout.runs[3].text=="x"&&layout.runs[3].y>0&&layout.runs[3].scale<1,
            "Actual font layout changed atom identity placement or lost AO script");
    }
    for(const std::string label:{"≈t₂g","2×t₂g"}) {
        const auto layout=cov::ui::scientific_text_layout(label,24);
        require(layout.runs.size()==3&&layout.runs[0].y==0&&layout.runs[0].scale==1&&
            layout.runs[1].text=="t"&&layout.runs[1].y==0&&layout.runs[2].text=="2g"&&
            layout.runs[2].y>0&&layout.runs[2].scale<1,"Prefixed irrep has incorrect run placement");
    }
    const auto approximate=cov::ui::scientific_text_layout("≈σg+",24);
    require(approximate.runs.size()==4&&approximate.runs[0].y==0&&approximate.runs[0].scale==1&&
        approximate.runs[2].text=="g"&&approximate.runs[2].y>0&&
        approximate.runs[3].text=="+"&&approximate.runs[3].y<0&&
        near(approximate.runs[2].x,approximate.runs[3].x),"Prefixed linear irrep lost its two script lanes");
    const auto linear=cov::ui::scientific_text_layout("2σg+",24);
    require(linear.runs.size()==3&&linear.runs[0].text=="2σ"&&
        linear.runs[1].text=="g"&&linear.runs[2].text=="+"&&
        near(linear.runs[1].x,linear.runs[2].x)&&linear.runs[1].y>0&&linear.runs[2].y<0,
        "Linear parity and reflection sign do not share opposite script lanes");
    const float expected=font->CalcTextSizeA(24,1e6f,0,"2σ").x+
        std::max(font->CalcTextSizeA(24*.72f,1e6f,0,"g").x,
                 font->CalcTextSizeA(24*.72f,1e6f,0,"+").x);
    require(near(expected,linear.width),"Script lanes incorrectly add both widths");
    const auto unicode=cov::ui::scientific_text_layout("t₂g",24);
    const auto ascii=cov::ui::scientific_text_layout("t2g",24);
    require(near(unicode.width,ascii.width)&&unicode.runs.size()==ascii.runs.size()&&
        near(unicode.runs[1].scale,ascii.runs[1].scale),"Unicode subscript was lowered or scaled twice");
    const auto spans=cov::scientific_svg_spans(linear,30,50,24);
    require(spans.find("y=\"58.64\"")!=std::string::npos&&
        spans.find("y=\"48.08\"")!=std::string::npos&&
        spans.find("font-size=\"17.28\"")!=std::string::npos,"SVG lost script coordinates or font size");
    require(spans.find("xml:space=\"preserve\"")!=std::string::npos&&
        spans.find("textLength=\"")!=std::string::npos&&
        spans.find("lengthAdjust=\"spacingAndGlyphs\"")!=std::string::npos,
        "SVG must preserve spaces and measured run widths across browser fonts");

    ImDrawList geometry(ImGui::GetDrawListSharedData());
    // Standalone lists need the same command-buffer initialization that ImGui
    // performs for window lists before clip or texture state can be pushed.
    geometry._ResetForNewFrame();geometry.PushClipRectFullScreen();
    geometry.PushTextureID(io.Fonts->TexID);
    cov::ui::draw_scientific_text(&geometry,24,{30,50},IM_COL32(20,30,40,255),"2σg+");
    require(geometry.VtxBuffer.Size==16,"Scientific text did not draw the expected glyphs");
    const auto base=bounds(geometry,0,8),sub=bounds(geometry,8,12),sup=bounds(geometry,12,16);
    require(sub.bottom>base.bottom&&sup.top<sub.top&&
        sub.left>=base.right-2&&sup.left>=base.right-2,
        "Actual ImGui glyph vertices do not occupy upper and lower script lanes");
    require(sub.bottom-sub.top<base.bottom-base.top,"Actual script glyph is not smaller than the base");

    ImDrawList image(ImGui::GetDrawListSharedData());
    image._ResetForNewFrame();image.PushClipRectFullScreen();
    image.PushTextureID(io.Fonts->TexID);
    const auto colour=IM_COL32(25,42,63,255);
    image.AddText(font,20,{24,14},colour,"Orbital and point-group typography (application font)");
    const std::string lines[]={
        "12t₂g [alpha]     3a₁′     E2''     2σg+     Σ_u^-     πᵤ",
        "C2v     D3h     Oh     D∞h     Cinfv     Cs     Ci     SO(3)",
        "O1 2px     Cr1 dxy / NAO 1     10p_x     d_x2-y2     dz2",
        "Gaussian AO 8 / shell 3 d x^2     f x^2y     g x^4",
        "MO 8 [beta] / NAO 12     d real m=1 cos     mixed directions"
    };
    for(std::size_t i=0;i<5;++i)
        cov::ui::draw_scientific_text(&image,28,{24,58+static_cast<float>(i)*57},colour,lines[i]);
    const auto path=argc>1?std::filesystem::path(argv[1]):std::filesystem::path("scientific_text_preview.ppm");
    preview(image,path,1100,365);
    auto bundle=path;bundle.replace_extension();bundle+="-bundle";
    export_typography(bundle);
    ImGui::EndFrame();ImGui::DestroyContext();
    std::cout<<"scientific_text_smoke passed; preview="<<path.string()<<'\n';return 0;
} catch(const std::exception& error) {
    std::cerr<<"scientific_text_smoke: "<<error.what()<<'\n';
    if(ImGui::GetCurrentContext())ImGui::DestroyContext();return 1;
} }
