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

cmake_minimum_required(VERSION 3.28)
project(cov_compute_sycl LANGUAGES CXX)

if(NOT COV_ROOT)
    get_filename_component(COV_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
endif()
if(NOT EXISTS "${COV_ROOT}/include/cov/compute_abi.h")
    message(FATAL_ERROR "COV_ROOT must point to the COV source root")
endif()
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    message(FATAL_ERROR "The SYCL module requires the Intel oneAPI DPC++ compiler (icpx)")
endif()

add_library(cov_compute_sycl SHARED compute_sycl.cpp)
target_include_directories(cov_compute_sycl PRIVATE "${COV_ROOT}/include")
target_compile_features(cov_compute_sycl PRIVATE cxx_std_17)
target_compile_options(cov_compute_sycl PRIVATE -fsycl)
target_link_options(cov_compute_sycl PRIVATE -fsycl)
set_target_properties(cov_compute_sycl PROPERTIES CXX_VISIBILITY_PRESET hidden)
