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

#pragma once

#include <algorithm>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace cov {

// Presentation only. These runs never replace source labels, orbital IDs or
// point-group identifiers. Both ASCII source names and Unicode display names
// are accepted, so a Unicode subscript is not scaled/lowered a second time.
struct ScientificTextPart {
    std::string base, subscript, superscript;
};
struct ScientificTextRun {
    std::string text;
    float x = 0, y = 0, scale = 1;
    float width = 0;
};
struct ScientificTextLayout {
    std::vector<ScientificTextRun> runs;
    float width = 0;
};

inline void scientific_replace(std::string& text, std::string_view from,
                               std::string_view to) {
    for (std::size_t at = 0; (at = text.find(from, at)) != std::string::npos;) {
        text.replace(at, from.size(), to);
        at += to.size();
    }
}

inline std::vector<ScientificTextPart> scientific_token(const std::string& raw,
                                                       bool point_group=false,
                                                       bool atom_indices=false) {
    // The caller knows whether B1 is an atom identity or an irrep. Preserve
    // literal atom IDs in atomic labels; explicitly scripted B₁ remains an irrep.
    static const std::regex atom_id(R"(^[A-Z][a-z]?[0-9]+$)");
    if(atom_indices && std::regex_match(raw,atom_id))return {{raw,{},{}}};
    // Approximation and copy-count prefixes are baseline text, not part of
    // the irrep suffix. Keep them when parsing the following scientific token.
    std::size_t prefix=raw.starts_with("≈")?std::string_view("≈").size():0;
    if(!prefix) {
        std::size_t digits=0;
        while(digits<raw.size() && raw[digits]>='0' && raw[digits]<='9')++digits;
        if(digits && std::string_view(raw).substr(digits).starts_with("×"))
            prefix=digits+std::string_view("×").size();
    }
    if(prefix && prefix<raw.size()) {
        auto parts=scientific_token(raw.substr(prefix),point_group,atom_indices);
        parts.insert(parts.begin(),{raw.substr(0,prefix),{},{}});
        return parts;
    }
    std::string token = raw;
    const std::pair<const char*, const char*> scripts[] = {
        {"₀","0"},{"₁","1"},{"₂","2"},{"₃","3"},{"₄","4"},
        {"₅","5"},{"₆","6"},{"₇","7"},{"₈","8"},{"₉","9"},
        {"ᵤ","u"},{"ₓ","x"},{"⁺","+"},{"⁻","-"},{"−","-"},
        {"′","'"},{"″","''"}
    };
    for (const auto& [from,to] : scripts) scientific_replace(token,from,to);
    std::smatch m;
    static const std::regex finite(R"(^([0-9]*[AaBbEeTtFfGgHh])_?([0-9]*[gGuU]?)(''|'|"|\*)?$)");
    if (std::regex_match(token,m,finite)) {
        std::string upper=m[3];
        scientific_replace(upper,"\"","′′");
        scientific_replace(upper,"'","′");
        return {{m[1],m[2],upper}};
    }
    // Keep the existing case: an orbital sigma and a state/irrep Sigma are
    // different uses of the same angular symbol.
    for (const auto symbol : {"σ","π","δ","φ","γ","Σ","Π","Δ","Φ","Γ"}) {
        std::size_t start=0;
        while(start<token.size() && token[start]>='0' && token[start]<='9')++start;
        if(!std::string_view(token).substr(start).starts_with(symbol))continue;
        const auto end=start+std::string_view(symbol).size();
        const std::string tail=token.substr(end);
        static const std::regex linear(R"(^_?([guGU]?)(?:\^?([+\-*]))?$)");
        if(std::regex_match(tail,m,linear)) {
            const std::string sign=m[2]=="-"?"−":m[2].str();
            return {{token.substr(0,end),m[1],sign}};
        }
    }
    static const std::regex group(R"(^([CDS])_?([0-9]+|n|inf|∞)([hvd]?)$)");
    if(std::regex_match(token,m,group) &&
       (point_group || !m[3].str().empty() || m[2]=="∞" || m[2]=="inf" || m[2]=="n"))
        return {{m[1],(m[2]=="inf"?std::string("∞"):m[2].str())+m[3].str(),{}}};
    static const std::regex special_group(R"(^([CTOI])_?([sihd])$)");
    if(std::regex_match(token,m,special_group))return {{m[1],m[2],{}}};

    // p_x, d_xy, d_x2-y2, etc. The principal quantum number is baseline.
    static const std::regex atomic(R"(^([0-9]*[spdfgh])_?([xyz][xyz0-9^+\-]*)$)");
    if(std::regex_match(token,m,atomic)) {
        std::string direction=m[2];
        direction.erase(std::remove(direction.begin(),direction.end(),'^'),direction.end());
        const char* powers[]={"⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹"};
        std::string sub;
        for(char c:direction) {
            if(c>='0' && c<='9')sub+=powers[c-'0'];
            else if(c=='-')sub+="−";
            else sub+=c;
        }
        return {{m[1],sub,{}}};
    }
    // Gaussian Cartesian polynomial powers, e.g. x^2y, have no AO
    // subscript wrapper. Leave spherical "real m=1 cos" and indices alone.
    static const std::regex polynomial(R"(^[xyz](?:\^[0-9]+)?(?:[xyz](?:\^[0-9]+)?)*$)");
    if(token.find('^')!=std::string::npos && std::regex_match(token,polynomial)) {
        std::vector<ScientificTextPart> result;
        for(std::size_t i=0;i<token.size();) {
            ScientificTextPart part;part.base=token[i++];
            if(i<token.size() && token[i]=='^') {
                ++i;
                while(i<token.size() && token[i]>='0' && token[i]<='9')
                    part.superscript+=token[i++];
            }
            result.push_back(std::move(part));
        }
        return result;
    }
    return {{raw,{},{}}};
}

inline std::vector<ScientificTextPart> scientific_text_parts(std::string_view text,
        bool point_group=false,bool atom_indices=false) {
    std::vector<ScientificTextPart> parts;
    const auto delimiter=[](char c) {
        return c==' ' || c=='\t' || c=='\n' || c==',' || c==';' ||
               c==':' || c=='/' || c=='[' || c==']' || c=='(' || c==')';
    };
    for(std::size_t start=0;start<text.size();) {
        std::size_t end=start+1;
        if(delimiter(text[start])) {
            parts.push_back({std::string(text.substr(start,1)),{},{}});
        } else {
            while(end<text.size() && !delimiter(text[end]))++end;
            auto token=scientific_token(std::string(text.substr(start,end-start)),point_group,atom_indices);
            parts.insert(parts.end(),token.begin(),token.end());
        }
        start=end;
    }
    return parts;
}

// measure(text, relative size) is supplied by the actual font renderer.
// Script lanes share a horizontal origin, as in sigma_g^+.
template<class Measure>
ScientificTextLayout layout_scientific_text(std::string_view text,float size,Measure measure,
        bool point_group=false,bool atom_indices=false) {
    ScientificTextLayout result;
    constexpr float script_scale=.72f;
    for(const auto& part:scientific_text_parts(text,point_group,atom_indices)) {
        if(!part.base.empty()) {
            const auto width=measure(part.base,1.0f);
            result.runs.push_back({part.base,result.width,0,1,width});
            result.width+=width;
        }
        float suffix=0;
        if(!part.subscript.empty()) {
            suffix=measure(part.subscript,script_scale);
            result.runs.push_back({part.subscript,result.width,.36f*size,script_scale,suffix});
        }
        if(!part.superscript.empty()) {
            const auto width=measure(part.superscript,script_scale);
            result.runs.push_back({part.superscript,result.width,-.08f*size,script_scale,width});
            suffix=std::max(suffix,width);
        }
        result.width+=suffix;
    }
    return result;
}

inline std::string scientific_xml(std::string_view value) {
    std::string out;
    for(char c:value) {
        switch(c) {
        case '&':out+="&amp;";break;case '<':out+="&lt;";break;
        case '>':out+="&gt;";break;case '"':out+="&quot;";break;
        default:out+=c;break;
        }
    }
    return out;
}
inline std::string scientific_svg_spans(const ScientificTextLayout& layout,
                                      float x,float y,float size) {
    std::ostringstream out;
    for(const auto& run:layout.runs) {
        out<<"<tspan x=\""<<x+run.x<<"\" y=\""<<y+run.y
           <<"\" font-size=\""<<size*run.scale<<"\" xml:space=\"preserve\"";
        // SVG may use a different installed font. Freeze each run's measured
        // width (including spaces), not the whole label with script lanes.
        if(run.width>0)out<<" textLength=\""<<run.width<<"\" lengthAdjust=\"spacingAndGlyphs\"";
        out<<">"
           <<scientific_xml(run.text)<<"</tspan>";
    }
    return out.str();
}
} // namespace cov
