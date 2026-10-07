if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    cmake_minimum_required(VERSION 3.21)
    project(cov_compute_webgpu_module LANGUAGES CXX)
endif()

if(NOT DEFINED COV_ROOT)
    get_filename_component(COV_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
endif()

# Pin the C ABI: other wgpu-native releases are not assumed source compatible.
# Windows x64 MSVC release archive:
# https://github.com/gfx-rs/wgpu-native/releases/download/v29.0.1.1/wgpu-windows-x86_64-msvc-release.zip
# SHA256 7e67d7445c42aeb85e30f88930fd8d7d83ee769e3390aeb1ada75ebf3cf78132
set(COV_WGPU_NATIVE_VERSION "29.0.1.1")
set(COV_WGPU_NATIVE_ROOT "" CACHE PATH "Extracted wgpu-native v${COV_WGPU_NATIVE_VERSION} package root")
if(NOT COV_WGPU_NATIVE_ROOT)
    message(FATAL_ERROR "COV_ENABLE_WEBGPU requires COV_WGPU_NATIVE_ROOT (wgpu-native v${COV_WGPU_NATIVE_VERSION})")
endif()

find_path(COV_WGPU_INCLUDE_DIR webgpu/wgpu.h PATHS "${COV_WGPU_NATIVE_ROOT}/include" NO_DEFAULT_PATH REQUIRED)
if(WIN32)
    find_library(COV_WGPU_LIBRARY NAMES wgpu_native.dll PATHS "${COV_WGPU_NATIVE_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
    find_file(COV_WGPU_RUNTIME NAMES wgpu_native.dll PATHS "${COV_WGPU_NATIVE_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
elseif(APPLE)
    find_file(COV_WGPU_RUNTIME NAMES libwgpu_native.dylib PATHS "${COV_WGPU_NATIVE_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
    set(COV_WGPU_LIBRARY "${COV_WGPU_RUNTIME}")
else()
    find_file(COV_WGPU_RUNTIME NAMES libwgpu_native.so PATHS "${COV_WGPU_NATIVE_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
    set(COV_WGPU_LIBRARY "${COV_WGPU_RUNTIME}")
endif()

file(READ "${CMAKE_CURRENT_LIST_DIR}/orbital.wgsl" COV_WGSL_SOURCE)
configure_file("${CMAKE_CURRENT_LIST_DIR}/shader.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/shader.hpp" @ONLY)

add_library(cov_compute_webgpu SHARED webgpu.cpp)
target_compile_features(cov_compute_webgpu PRIVATE cxx_std_17)
target_include_directories(cov_compute_webgpu PRIVATE "${COV_ROOT}/include" "${COV_WGPU_INCLUDE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(cov_compute_webgpu PRIVATE "${COV_WGPU_LIBRARY}")
add_custom_command(TARGET cov_compute_webgpu POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${COV_WGPU_RUNTIME}" "$<TARGET_FILE_DIR:cov_compute_webgpu>")
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(COV_WGPU_RUNTIME "${COV_WGPU_RUNTIME}" PARENT_SCOPE)
endif()
if(APPLE)
    set_target_properties(cov_compute_webgpu PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "@loader_path")
elseif(UNIX)
    set_target_properties(cov_compute_webgpu PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH "$ORIGIN")
endif()
