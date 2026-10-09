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
