cmake_minimum_required(VERSION 3.28)

# The Windows HIP SDK does not implement CMake's HIP language. hipcc compiles
# HIP sources as C++ there; Linux uses CMake's native HIP support.
if(CMAKE_HOST_WIN32)
    project(cov_compute_hip LANGUAGES CXX)
    if(NOT CMAKE_CXX_COMPILER MATCHES "(^|[/\\\\])hipcc(\\.exe)?$")
        message(FATAL_ERROR "On Windows configure this module with -DCMAKE_CXX_COMPILER=<HIP SDK>/bin/hipcc.exe")
    endif()
else()
    project(cov_compute_hip LANGUAGES HIP)
    set_source_files_properties(compute_hip.cpp PROPERTIES LANGUAGE HIP)
endif()

if(NOT COV_ROOT)
    get_filename_component(COV_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
endif()
if(NOT EXISTS "${COV_ROOT}/include/cov/compute_abi.h")
    message(FATAL_ERROR "COV_ROOT must point to the COV source root")
endif()

add_library(cov_compute_hip SHARED compute_hip.cpp)
target_include_directories(cov_compute_hip PRIVATE "${COV_ROOT}/include")
set_target_properties(cov_compute_hip PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    HIP_VISIBILITY_PRESET hidden
)
if(CMAKE_HOST_WIN32)
    target_compile_features(cov_compute_hip PRIVATE cxx_std_17)
    if(NOT COV_HIP_ARCHITECTURES)
        message(FATAL_ERROR "Set COV_HIP_ARCHITECTURES to the target gfx architecture on Windows")
    endif()
    foreach(arch IN LISTS COV_HIP_ARCHITECTURES)
        target_compile_options(cov_compute_hip PRIVATE "--offload-arch=${arch}")
        target_link_options(cov_compute_hip PRIVATE "--offload-arch=${arch}")
    endforeach()
else()
    set_target_properties(cov_compute_hip PROPERTIES HIP_STANDARD 17 HIP_STANDARD_REQUIRED ON)
endif()
