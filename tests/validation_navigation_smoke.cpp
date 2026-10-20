#include "cov/validation_navigation.hpp"
#include <imgui_internal.h>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition,const char* message) { if (!condition) throw std::runtime_error(message); }
void nested_navigation(float scale,bool reveal_entire_item=false,bool fractional_vertical=false) {
    ImGui::CreateContext();auto& io=ImGui::GetIO();
    io.IniFilename=nullptr;io.DisplaySize={1200,1100};io.DeltaTime=1.0f/60;
    io.Fonts->AddFontDefault();io.Fonts->Build();
    io.ConfigInputTrickleEventQueue=true;
    std::array<cov::validation::NavigationTarget,4> targets;
    ImGuiWindow* canvas=nullptr;ImGuiWindow* card=nullptr;
    int selected=-1;ImVec2 pointer(-100,-100);
    auto draw=[&]() {
        ImGui::NewFrame();ImGui::SetNextWindowPos({10,10},ImGuiCond_Always);
        ImGui::SetNextWindowSize({460*scale,620*scale},ImGuiCond_Always);
        ImGui::Begin("navigation root",nullptr,ImGuiWindowFlags_NoSavedSettings);
        ImGui::Dummy({0,30*scale});
        ImGui::BeginChild("outer panel",{0,0});
        ImGui::Dummy({0,110*scale});
        ImGui::BeginChild("card",{0,310*scale},true);
        card=ImGui::GetCurrentWindow();ImGui::Dummy({0,45*scale});
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding,{0,0});
        ImGui::BeginChild("canvas",{0,430*scale+ImGui::GetStyle().ScrollbarSize},
                          ImGuiChildFlags_None,ImGuiWindowFlags_HorizontalScrollbar);
        ImGui::PopStyleVar();canvas=ImGui::GetCurrentWindow();
        const auto origin=ImGui::GetCursorScreenPos();
        ImGui::Dummy({1200*scale,(fractional_vertical?1200:430)*scale});
        const float top=fractional_vertical?85.5f:85.0f,bottom=fractional_vertical?955.5f:355.0f;
        const std::array<ImVec2,4> offsets{{{950,top},{60,top},{950,bottom},{60,bottom}}};
        for (int i=0;i<4;++i) {
            ImGui::SetCursorScreenPos({origin.x+offsets[i].x*scale,origin.y+offsets[i].y*scale});
            ImGui::PushID(i);
            if (ImGui::Button("target",{50*scale,20*scale})) selected=i;
            targets[i]={ImGui::GetItemRectMin(),ImGui::GetItemRectMax(),canvas,canvas->ClipRect.Min,canvas->ClipRect.Max};
            ImGui::PopID();
        }
        ImGui::EndChild();ImGui::EndChild();
        ImGui::Dummy({0,500*scale});ImGui::EndChild();ImGui::End();ImGui::Render();
    };
    auto tick=[&](ImVec2 wheel=ImVec2(0,0),int button=-1) {
        // Same per-frame cleanup and ordinary input sequence as native_session.
        io.ClearEventsQueue();io.AddFocusEvent(true);io.AddMousePosEvent(pointer.x,pointer.y);
        if (wheel.x || wheel.y) io.AddMouseWheelEvent(wheel.x,wheel.y);
        if (button>=0) io.AddMouseButtonEvent(0,button!=0);
        draw();
    };
    tick();tick();tick();
    require(canvas->ScrollMax.x>0,"fixture must need horizontal scrolling");
    const auto bar=ImGui::GetWindowScrollbarRect(canvas,ImGuiAxis_X);
    require(bar.Min.y>card->ClipRect.Max.y,"fixture must clip the inner horizontal scrollbar");
    require(targets[0].lo.x>targets[0].clip_hi.x,"first target must start outside the horizontal clip");
    require(targets[0].lo.y>=targets[0].clip_lo.y && targets[0].hi.y<=targets[0].clip_hi.y,
            "first target must already be vertically visible");
    int moves=0,wheels_x=0,wheels_y=0,visit=0;
    for (int index:{0,1,2,3,0}) {
        bool reached=false;
        for (int attempt=0;attempt<60;++attempt) {
            const auto step=cov::validation::plan_navigation(targets[index],reveal_entire_item);
            if (attempt==59) {
                std::cerr<<"scale="<<scale<<" target="<<index<<" step="<<int(step.kind)
                    <<" mouse="<<io.MousePos.x<<','<<io.MousePos.y
                    <<" requested="<<step.mouse.x<<','<<step.mouse.y
                    <<" canvas_scroll="<<canvas->Scroll.x<<','<<canvas->Scroll.y
                    <<" card_scroll="<<card->Scroll.x<<','<<card->Scroll.y<<'\n';
            }
            require(step.kind!=cov::validation::NavigationKind::Unreachable,"nested target unreachable through actual wheel input");
            pointer=step.mouse;
            if (step.kind==cov::validation::NavigationKind::Ready) {
                if(reveal_entire_item) {
                    const auto& t=targets[index];
                    require(t.lo.x>=t.clip_lo.x && t.lo.y>=t.clip_lo.y &&
                            t.hi.x<=t.clip_hi.x && t.hi.y<=t.clip_hi.y,
                            "hover evidence must reveal the complete fitting target in the actual clip");
                }
                if (visit==0) require(wheels_y==0,"a horizontal miss must not realign vertical ancestors");
                tick();tick({},1);tick({},0);
                require(selected==index,"actual button did not receive the selected target click");
                reached=true;break;
            }
            if (step.kind==cov::validation::NavigationKind::Move) { ++moves;tick(); }
            else {
                require(std::abs(io.MousePos.x-pointer.x)<0.01f && std::abs(io.MousePos.y-pointer.y)<0.01f,
                        "wheel and pointer move must not share the trickled frame");
                wheels_x+=step.wheel.x!=0;wheels_y+=step.wheel.y!=0;
                tick(step.wheel);tick();tick();tick();
            }
        }
        require(reached,"navigation exceeded its existing bounded attempt count");
        ++visit;
    }
    require(moves>0 && wheels_x>0 && wheels_y>0,"fixture must exercise both axes and staged pointer moves");
    ImGui::DestroyContext();
}
}
int main() {
    try {
        nested_navigation(1);nested_navigation(1.5f);
        nested_navigation(1,true);nested_navigation(1.5f,true);
        nested_navigation(1,true,true);nested_navigation(1.5f,true,true);
        std::cout<<"validation_navigation_smoke ok\n";
    } catch (const std::exception& error) {
        std::cerr<<"validation_navigation_smoke: "<<error.what()<<'\n';
        return 1;
    }
}
