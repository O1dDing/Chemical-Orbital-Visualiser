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

#include "cov/viewer_layout.hpp"
#include <iomanip>
#include <iostream>

int main() {
    std::cout << std::setprecision(17);
    int id = 0, fw = 0, fh = 0;
    float w = 0, h = 0, scale = 0, fov = 0;
    while (std::cin >> id >> w >> h >> fw >> fh >> scale >> fov) {
        const auto layout = cov::viewer_layout(w, h, fw, fh, scale);
        const auto& s = layout.scene;
        const auto& c = layout.controls;
        const auto& v = layout.framebuffer;
        const auto p = cov::scene_projection(v.width, v.height, fov);
        std::cout << "{\"id\":" << id << ",\"controls\":[" << c.x << ',' << c.y << ',' << c.width << ',' << c.height
            << "],\"scene\":[" << s.x << ',' << s.y << ',' << s.width << ',' << s.height
            << "],\"viewport\":[" << v.x << ',' << v.y << ',' << v.width << ',' << v.height
            << "],\"projection\":[" << p.aspect << ',' << p.tan_half_vertical_fov
            << "],\"hit\":[" << s.contains(s.x+s.width*.5f,s.y+s.height*.5f) << ','
            << s.contains(s.x+s.width,s.y+s.height*.5f) << ','
            << s.contains(c.x+c.width*.5f,c.y+c.height*.5f) << "]}\n";
    }
    return std::cin.eof() ? 0 : 2;
}
