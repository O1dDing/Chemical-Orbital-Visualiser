# SPDX-License-Identifier: CPAL-1.0
#
# The contents of this file are subject to the Common Public Attribution License
# Version 1.0 (the “License”); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
# https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
# License Version 1.1 but Sections 14 and 15 have been added to cover use of software
# over a computer network and provide for limited attribution for the Original
# Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
#
# Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
# WARRANTY OF ANY KIND, either express or implied. See the License for the specific
# language governing rights and limitations under the License.
#
# The Original Code is Chemical Orbital Visualiser (COV).
#
# The Original Developer is O1dDing, who is also the Initial Developer.
#
# The Initial Developer of the Original Code is O1dDing. All portions of the code
# written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
#
# Contributor O1dDing.

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    cmake_minimum_required(VERSION 3.21)
    project(cov_compute_metal_module LANGUAGES CXX OBJCXX)
endif()
if(NOT APPLE)
    message(FATAL_ERROR "Metal compute is available only on Apple platforms")
endif()
enable_language(OBJCXX)
if(NOT DEFINED COV_ROOT)
    get_filename_component(COV_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
endif()

set(COV_METAL_AIR "${CMAKE_CURRENT_BINARY_DIR}/orbital.air")
set(COV_METAL_LIB "${CMAKE_CURRENT_BINARY_DIR}/orbital.metallib")
set(COV_METAL_DEPLOYMENT_FLAGS)
if(CMAKE_OSX_DEPLOYMENT_TARGET)
    list(APPEND COV_METAL_DEPLOYMENT_FLAGS "-mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
endif()
add_custom_command(OUTPUT "${COV_METAL_LIB}"
    COMMAND xcrun -sdk macosx metal ${COV_METAL_DEPLOYMENT_FLAGS} -c "${CMAKE_CURRENT_LIST_DIR}/orbital.metal" -o "${COV_METAL_AIR}"
    COMMAND xcrun -sdk macosx metallib "${COV_METAL_AIR}" -o "${COV_METAL_LIB}"
    DEPENDS "${CMAKE_CURRENT_LIST_DIR}/orbital.metal"
    VERBATIM)
add_custom_target(cov_metal_shader DEPENDS "${COV_METAL_LIB}")

add_library(cov_compute_metal SHARED metal.mm)
add_dependencies(cov_compute_metal cov_metal_shader)
target_compile_features(cov_compute_metal PRIVATE cxx_std_17)
target_compile_options(cov_compute_metal PRIVATE -fobjc-arc)
target_include_directories(cov_compute_metal PRIVATE "${COV_ROOT}/include")
find_library(COV_METAL_FRAMEWORK Metal REQUIRED)
find_library(COV_FOUNDATION_FRAMEWORK Foundation REQUIRED)
target_link_libraries(cov_compute_metal PRIVATE "${COV_METAL_FRAMEWORK}" "${COV_FOUNDATION_FRAMEWORK}")
# The metallib is a Mach-O section in this dylib, not a file looked up at runtime.
target_link_options(cov_compute_metal PRIVATE "-Wl,-sectcreate,__TEXT,__cov_metal,${COV_METAL_LIB}")
