cmake_minimum_required(VERSION 3.28)

project(CUDAOrbitalVisualisation VERSION 0.4.0 LANGUAGES CXX)

option(COV_ENABLE_CUDA "Build the CUDA/OpenGL viewer" ON)
option(COV_BUILD_TESTS "Build parser, orbital, diagram, molecule and UI smoke tests" ON)
option(COV_ENABLE_VALIDATION "Enable local native validation plan and evidence capture" OFF)
if(COV_ENABLE_VALIDATION AND NOT COV_ENABLE_CUDA)
    message(FATAL_ERROR "Native validation requires the production CUDA/OpenGL viewer")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Keep the downloadable Windows viewer self-contained with respect to the
# Microsoft C/C++ runtime. OpenGL and the NVIDIA display driver remain normal
# operating-system/driver dependencies.
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

add_library(cov_parser
    src/parser/molden_parser.cpp
    src/parser/molden_compat.cpp
    src/parser/fchk_parser.cpp
    src/parser/fchk_overlap.cpp
    src/parser/fchk_nuclear_charge.cpp
    src/parser/wavefunction_io.cpp
    src/parser/formchk.cpp
    src/parser/gaussian_log.cpp
    src/parser/nbo.cpp
    src/parser/nbo_integration.cpp
    src/analysis/nbo_salc.cpp
    src/analysis/nbo_spin_average.cpp
    src/analysis/pi_coupling.cpp
    src/analysis/pi_field_response.cpp
    src/analysis/density.cpp
    src/analysis/overlap.cpp
    src/analysis/overlap_basis.cpp
    src/metric_subspace.cpp
    src/local_angular_generators.cpp
    src/analysis/atomic_orbital_symmetry.cpp
    src/local_angular_projection.cpp
    src/analysis/bond_analysis.cpp
    src/analysis/symmetry.cpp
    src/analysis/point_group_catalog.cpp
    src/analysis/point_group_irreps.cpp
    src/analysis/orbital_symmetry.cpp
    src/analysis/orbital_symmetry_final.cpp
    src/analysis/orbital_chemistry.cpp
    src/analysis/pi_topology.cpp
    src/analysis/topology_graph.cpp
)
set_source_files_properties(src/parser/molden_parser.cpp PROPERTIES
    COMPILE_DEFINITIONS "parse_molden=parse_molden_format"
)
set_source_files_properties(src/analysis/orbital_symmetry.cpp PROPERTIES
    COMPILE_DEFINITIONS "derive_orbital_symmetry=derive_orbital_symmetry_legacy"
)
target_include_directories(cov_parser PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(cov_parser PUBLIC cxx_std_20)

# Header-only dense spectral diagnostics. Pin the source archive and checksum;
# do not enable Eigen's tests or add a runtime DLL dependency to cov.exe.
include(FetchContent)
FetchContent_Declare(cov_eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR cov_headers_only
)
FetchContent_MakeAvailable(cov_eigen)
target_include_directories(cov_parser SYSTEM PRIVATE ${cov_eigen_SOURCE_DIR})
target_compile_definitions(cov_parser PRIVATE EIGEN_MPL2_ONLY)

add_library(cov_orbitals
    src/analysis/nbo_channels.cpp
    src/analysis/local_orbital_symmetry.cpp
    src/analysis/local_geometry.cpp
    src/analysis/orbital_tracking.cpp
    src/orbitals/orbital_view.cpp
    src/orbitals/orbital_view_dispatch.cpp
    src/orbitals/orbital_symmetry_scope.cpp
    src/orbitals/orbital_symmetry_json.cpp
    src/orbitals/pi_pair_evidence.cpp
    src/orbitals/orbital_energy_gap.cpp
    src/orbitals/coordination_geometry.cpp
    src/orbitals/ligand_field.cpp
    src/orbitals/mo_diagram_v2.cpp
    src/orbitals/mo_diagram_snapshot.cpp
    src/orbitals/mo_diagram_layout.cpp
    src/orbitals/mo_diagram_chemistry.cpp
    src/orbitals/mo_sigma_framework.cpp
    src/orbitals/annotation_source_dispatch.cpp
)
set_source_files_properties(src/orbitals/orbital_view.cpp PROPERTIES
    COMPILE_DEFINITIONS "classify_orbital_region=classify_orbital_region_legacy;orbital_visible=orbital_visible_legacy;visible_orbital_indices=visible_orbital_indices_legacy;build_orbital_metadata=build_orbital_metadata_legacy"
)
set_source_files_properties(src/orbitals/mo_diagram_v2.cpp PROPERTIES
    COMPILE_DEFINITIONS "annotation_source_name=annotation_source_name_legacy;annotate_orbital=annotate_orbital_legacy;build_valence_selection_plan=build_valence_selection_plan_legacy;build_mo_diagram_data=build_mo_diagram_data_legacy"
)
target_include_directories(cov_orbitals PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(cov_orbitals PUBLIC cov_parser)
target_compile_features(cov_orbitals PUBLIC cxx_std_20)

add_library(cov_molecule
    src/analysis/orbital_group_bonding.cpp
    src/analysis/chemistry_route.cpp
    src/analysis/interaction_graph.cpp
    src/render/molecule_style.cpp
    src/render/nbo_molecular_overlay.cpp
)
target_include_directories(cov_molecule PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(cov_molecule PUBLIC cxx_std_20)
# Chemistry routing consumes parser-side NAO/SALC/Fock evidence. Express the
# dependency so static-library consumers link in the correct order on GNU ld.
target_link_libraries(cov_molecule PUBLIC cov_parser)
target_link_libraries(cov_orbitals PUBLIC cov_molecule)
add_library(cov_platform
    src/platform/file_dialog.cpp
)
target_include_directories(cov_platform PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_features(cov_platform PUBLIC cxx_std_20)
if(WIN32)
    target_link_libraries(cov_platform PRIVATE comdlg32)
endif()

if(MSVC)
    target_compile_options(cov_parser PRIVATE /W4 /permissive- /utf-8)
    target_compile_options(cov_orbitals PRIVATE /W4 /permissive- /utf-8)
    target_compile_options(cov_molecule PRIVATE /W4 /permissive- /utf-8)
    target_compile_options(cov_platform PRIVATE /W4 /permissive- /utf-8)
else()
    target_compile_options(cov_parser PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cov_orbitals PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cov_molecule PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(cov_platform PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(COV_ENABLE_CUDA OR COV_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        imgui
        URL https://codeload.github.com/ocornut/imgui/zip/refs/tags/v1.90.9
        URL_HASH SHA256=4fa247c7abbb00fc255631479004de549b86a41f49fddd7d44cd2cf8d647c5a7
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(imgui)

    add_library(cov_imgui_core STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    )
    target_include_directories(cov_imgui_core PUBLIC
        ${imgui_SOURCE_DIR}
    )
    # All consumers and backends must share the same index ABI. Large complete
    # AO/MO diagrams exceed 65535 vertices in a single ImGui draw list.
    target_compile_definitions(cov_imgui_core PUBLIC "ImDrawIdx=unsigned int")

    add_library(cov_validation_navigation STATIC src/validation/navigation.cpp)
    target_include_directories(cov_validation_navigation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(cov_validation_navigation PUBLIC cov_imgui_core)
    target_compile_features(cov_validation_navigation PUBLIC cxx_std_20)

    add_library(cov_ui STATIC
        src/ui/ui.cpp
        src/ui/nbo_ui.cpp
        src/ui/nbo_aomo_ui.cpp
        src/ui/nbo_aomo_labels.cpp
        src/ui/orbital_symmetry_components.cpp
        src/ui/nbo_aomo_hover.cpp
        src/ui/nbo_aomo_text.cpp
        src/ui/ui_text_dispatch.cpp
        src/ui/orbital_ui_text.cpp
        src/ui/orbital_ui_v2.cpp
        src/ui/orbital_ui_dispatch.cpp
    )
    set_source_files_properties(src/ui/ui.cpp PROPERTIES
        COMPILE_DEFINITIONS "tr=tr_legacy"
    )
    set_source_files_properties(src/ui/orbital_ui_v2.cpp PROPERTIES
        COMPILE_DEFINITIONS "draw_orbital_browser=draw_orbital_browser_legacy;draw_energy_diagram=draw_energy_diagram_legacy"
    )
    target_include_directories(cov_ui PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${imgui_SOURCE_DIR}
    )
    target_link_libraries(cov_ui PUBLIC cov_imgui_core cov_orbitals cov_molecule)
    target_compile_features(cov_ui PUBLIC cxx_std_20)
    target_include_directories(cov_ui SYSTEM PRIVATE ${cov_eigen_SOURCE_DIR})
    target_compile_options(cov_ui PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
endif()

if(COV_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR "COV_ENABLE_CUDA=ON but no CUDA compiler was found.")
    endif()

    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "120-real;120-virtual" CACHE STRING
            "CUDA architectures (default: Blackwell sm_120 + PTX)" FORCE)
    endif()

    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    find_package(CUDAToolkit 12.8 REQUIRED)
    find_package(OpenGL REQUIRED)

    # AO-overlap recovery is ordinary C++ code with an optional dense CUDA
    # backend. Link the solver/BLAS runtime only in CUDA builds; CPU-only
    # builds retain the same public API and multithreaded fallback.
    target_compile_definitions(cov_parser PRIVATE COV_ENABLE_CUDA=1)
    target_link_libraries(cov_parser PRIVATE CUDA::cudart_static)
    if(WIN32)
        # NVIDIA's Windows toolkit ships cuBLAS/cuSOLVER import libraries and
        # DLLs, not static archives. Resolve them optionally at runtime so the
        # single-file cov.exe still starts on a driver-only installation and
        # falls back to the CPU implementation when the toolkit DLLs are absent.
        target_compile_definitions(cov_parser PRIVATE
            COV_CUDA_DYNAMIC_DENSE=1)
    else()
        target_link_libraries(cov_parser PRIVATE
            CUDA::cublas_static
            CUDA::cusolver_static
        )
    endif()

    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        glfw
        URL https://codeload.github.com/glfw/glfw/zip/refs/tags/3.4
        URL_HASH SHA256=a133ddc3d3c66143eba9035621db8e0bcf34dba1ee9514a9e23e96afd39fd57a
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(glfw)

    add_library(cov_imgui_backends STATIC
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp
    )
    target_include_directories(cov_imgui_backends PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
    )
    target_link_libraries(cov_imgui_backends PUBLIC
        cov_imgui_core
        glfw
        OpenGL::GL
    )

    add_library(cov_gpu
        src/cuda/orbital.cu
        src/render/gl_api.cpp
        src/render/volume_renderer.cpp
    )
    target_include_directories(cov_gpu PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(cov_gpu
        PUBLIC
            cov_parser
            cov_molecule
            CUDA::cudart_static
            OpenGL::GL
            glfw
    )
    set_target_properties(cov_gpu PROPERTIES
        CUDA_SEPARABLE_COMPILATION OFF
        CUDA_RESOLVE_DEVICE_SYMBOLS OFF
        POSITION_INDEPENDENT_CODE ON
    )
    target_compile_options(cov_gpu PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>
    )

    add_executable(cov
        src/main.cpp
    )
    target_include_directories(cov PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(cov PRIVATE
        cov_gpu
        cov_orbitals
        cov_molecule
        cov_platform
        cov_ui
        cov_imgui_backends
    )
    target_compile_options(cov PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
    if(WIN32)
        set_property(TARGET cov PROPERTY WIN32_EXECUTABLE TRUE)
        if(MSVC)
            target_sources(cov PRIVATE src/platform/cov.manifest)
            target_link_options(cov PRIVATE /ENTRY:mainCRTStartup)
        else()
            target_link_options(cov PRIVATE -Wl,-e,mainCRTStartup)
        endif()
    endif()
    if(COV_ENABLE_VALIDATION)
        execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            OUTPUT_VARIABLE COV_VALIDATION_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND git diff --quiet WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE COV_VALIDATION_DIRTY)
        if(COV_VALIDATION_DIRTY)
            string(APPEND COV_VALIDATION_COMMIT "-dirty")
        endif()
        add_library(cov_validation_support STATIC src/validation/native_session.cpp src/validation/forensic_capture.cpp)
        target_include_directories(cov_validation_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_compile_definitions(cov_validation_support PUBLIC COV_ENABLE_VALIDATION=1)
        target_compile_definitions(cov_validation_support PRIVATE COV_VALIDATION_COMMIT="${COV_VALIDATION_COMMIT}")
        target_compile_options(cov_validation_support PRIVATE $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
        target_link_libraries(cov_validation_support PUBLIC cov_imgui_core cov_validation_navigation OpenGL::GL)
        target_link_libraries(cov_ui PUBLIC cov_validation_support)
        set_target_properties(cov PROPERTIES OUTPUT_NAME cov_validation)
        if(COV_BUILD_TESTS)
            add_executable(cov_forensic_capture_smoke tests/forensic_capture_smoke.cpp)
            target_link_libraries(cov_forensic_capture_smoke PRIVATE cov_validation_support)
            target_compile_options(cov_forensic_capture_smoke PRIVATE $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
            add_test(NAME forensic_capture_smoke COMMAND cov_forensic_capture_smoke)
        endif()
    endif()
endif()

if(COV_BUILD_TESTS)
    add_executable(cov_nbo_skeleton_policy_smoke tests/nbo_skeleton_policy_smoke.cpp)
    add_executable(cov_orbital_group_bonding_smoke tests/orbital_group_bonding_smoke.cpp)
    target_link_libraries(cov_orbital_group_bonding_smoke PRIVATE cov_molecule cov_orbitals cov_parser)
    add_test(NAME orbital_group_bonding_smoke COMMAND cov_orbital_group_bonding_smoke)
    add_executable(cov_orbital_group_bonding_probe tests/orbital_group_bonding_probe.cpp)
    target_link_libraries(cov_orbital_group_bonding_probe PRIVATE cov_molecule cov_orbitals cov_parser)
    target_link_libraries(cov_nbo_skeleton_policy_smoke PRIVATE cov_molecule cov_orbitals cov_parser)
    add_test(NAME nbo_skeleton_policy_smoke COMMAND cov_nbo_skeleton_policy_smoke)
    add_executable(cov_nbo_skeleton_policy_probe tests/nbo_skeleton_policy_probe.cpp)
    target_link_libraries(cov_nbo_skeleton_policy_probe PRIVATE cov_molecule cov_orbitals cov_parser)
    add_executable(cov_chemistry_route_smoke tests/chemistry_route_smoke.cpp)
    target_link_libraries(cov_chemistry_route_smoke PRIVATE cov_molecule cov_orbitals cov_parser)
    add_test(NAME chemistry_route_smoke COMMAND cov_chemistry_route_smoke)
    add_executable(cov_route_relation_projection_smoke tests/route_relation_projection_smoke.cpp)
    target_link_libraries(cov_route_relation_projection_smoke PRIVATE cov_molecule cov_orbitals cov_parser)
    add_test(NAME route_relation_projection_smoke COMMAND cov_route_relation_projection_smoke)
    add_executable(cov_route_relation_projection_probe tests/route_relation_projection_probe.cpp)
    target_link_libraries(cov_route_relation_projection_probe PRIVATE cov_molecule cov_orbitals cov_parser)
    add_executable(cov_chemistry_route_probe tests/chemistry_route_probe.cpp)
    target_link_libraries(cov_chemistry_route_probe PRIVATE cov_molecule cov_orbitals cov_parser)
    add_executable(cov_unified_aomo_probe tests/unified_aomo_probe.cpp)
    target_link_libraries(cov_unified_aomo_probe PRIVATE cov_parser)
    add_executable(cov_nbo_matrix_equivalence tests/nbo_matrix_equivalence.cpp)
    target_link_libraries(cov_nbo_matrix_equivalence PRIVATE cov_parser)
    add_executable(cov_nbo_salc_smoke tests/test_nbo_salc.cpp)
    add_executable(cov_nbo_spin_average_smoke tests/nbo_spin_average_smoke.cpp)
    target_link_libraries(cov_nbo_spin_average_smoke PRIVATE cov_parser)
    target_include_directories(cov_nbo_spin_average_smoke SYSTEM PRIVATE ${cov_eigen_SOURCE_DIR})
    add_test(NAME nbo_spin_average_smoke COMMAND cov_nbo_spin_average_smoke)
    add_executable(cov_nbo_spin_average_probe tests/nbo_spin_average_probe.cpp)
    target_link_libraries(cov_nbo_spin_average_probe PRIVATE cov_parser)
    add_executable(cov_pi_field_response_smoke tests/pi_field_response_smoke.cpp)
    target_link_libraries(cov_pi_field_response_smoke PRIVATE cov_parser)
    target_include_directories(cov_pi_field_response_smoke SYSTEM PRIVATE ${cov_eigen_SOURCE_DIR})
    add_test(NAME pi_field_response_smoke COMMAND cov_pi_field_response_smoke)
    add_executable(cov_pi_field_response_probe tests/pi_field_response_probe.cpp)
    target_link_libraries(cov_pi_field_response_probe PRIVATE cov_parser cov_orbitals cov_molecule)
    add_executable(cov_composition_display_policy_smoke tests/composition_display_policy_smoke.cpp)
    target_link_libraries(cov_composition_display_policy_smoke PRIVATE cov_parser cov_orbitals)
    add_test(NAME composition_display_policy_smoke COMMAND cov_composition_display_policy_smoke)
    target_link_libraries(cov_nbo_salc_smoke PRIVATE cov_parser)
    add_test(NAME nbo_salc_smoke COMMAND cov_nbo_salc_smoke)
    add_executable(cov_nbo_pi_evidence_smoke tests/test_nbo_pi_evidence.cpp)
    target_link_libraries(cov_nbo_pi_evidence_smoke PRIVATE cov_parser cov_orbitals)
    add_test(NAME nbo_pi_evidence_smoke COMMAND cov_nbo_pi_evidence_smoke)
    add_executable(cov_nbo_salc_invariance tests/nbo_salc_invariance.cpp)
    add_executable(cov_nbo_salc_copy_probe tests/nbo_salc_copy_probe.cpp)
    target_link_libraries(cov_nbo_salc_copy_probe PRIVATE cov_parser)
    target_include_directories(cov_nbo_salc_copy_probe SYSTEM PRIVATE ${cov_eigen_SOURCE_DIR})
    target_link_libraries(cov_nbo_salc_invariance PRIVATE cov_parser)
    target_include_directories(cov_nbo_salc_invariance SYSTEM PRIVATE ${cov_eigen_SOURCE_DIR})
    add_executable(cov_aomo_core_probe tests/aomo_core_probe.cpp)
    target_link_libraries(cov_aomo_core_probe PRIVATE cov_parser)
    add_executable(cov_nbo_integration_core_smoke tests/nbo_integration_core_smoke.cpp)
    target_link_libraries(cov_nbo_integration_core_smoke PRIVATE cov_parser cov_orbitals)
    add_test(NAME nbo_integration_core_smoke COMMAND cov_nbo_integration_core_smoke)
    include(CTest)
    enable_testing()

    add_executable(cov_nbo_probe tests/nbo_probe.cpp)
    target_link_libraries(cov_nbo_probe PRIVATE cov_parser)
    add_executable(cov_nbo_smoke tests/nbo_smoke.cpp)
    target_link_libraries(cov_nbo_smoke PRIVATE cov_parser)
    add_test(NAME nbo_smoke COMMAND cov_nbo_smoke)
    add_executable(cov_nbo_focus_core_smoke tests/nbo_focus_core_smoke.cpp)
    target_link_libraries(cov_nbo_focus_core_smoke PRIVATE cov_parser)
    add_test(NAME nbo_focus_core_smoke COMMAND cov_nbo_focus_core_smoke)

    if(COV_ENABLE_CUDA)
        # External scientific reference; never counted as a passing test when
        # the producer cube is absent. Exercises the production CUDA texture.
        add_executable(cov_cuda_cube_reference tests/cuda_cube_reference.cpp)
        target_link_libraries(cov_cuda_cube_reference PRIVATE cov_gpu)
    endif()

    add_executable(cov_parser_smoke tests/parser_smoke.cpp)
    target_link_libraries(cov_parser_smoke PRIVATE cov_parser)
    add_test(NAME parser_smoke
        COMMAND cov_parser_smoke ${CMAKE_CURRENT_SOURCE_DIR}/examples/h2.molden)

    add_executable(cov_fchk_parser_smoke tests/fchk_parser_smoke.cpp)
    target_link_libraries(cov_fchk_parser_smoke PRIVATE cov_parser)
    add_test(NAME fchk_parser_smoke COMMAND cov_fchk_parser_smoke)
    add_executable(cov_fchk_ao_contract_smoke tests/fchk_ao_contract_smoke.cpp)
    target_link_libraries(cov_fchk_ao_contract_smoke PRIVATE cov_parser)
    add_test(NAME fchk_ao_contract_smoke COMMAND cov_fchk_ao_contract_smoke)

    add_executable(cov_fchk_overlap_smoke tests/fchk_overlap_smoke.cpp)
    target_link_libraries(cov_fchk_overlap_smoke PRIVATE cov_parser)
    add_test(NAME fchk_overlap_smoke COMMAND cov_fchk_overlap_smoke)
    add_executable(cov_overlap_basis_smoke tests/overlap_basis_smoke.cpp)
    target_link_libraries(cov_overlap_basis_smoke PRIVATE cov_parser)
    add_test(NAME overlap_basis_smoke COMMAND cov_overlap_basis_smoke)

    add_executable(cov_format_equivalence_smoke tests/format_equivalence_smoke.cpp)
    target_link_libraries(cov_format_equivalence_smoke PRIVATE cov_parser)
    add_test(NAME format_equivalence_smoke COMMAND cov_format_equivalence_smoke)

    add_executable(cov_fchk_real_regression tests/fchk_real_regression.cpp)
    target_link_libraries(cov_fchk_real_regression PRIVATE cov_parser)
    add_test(NAME fchk_real_regression
        COMMAND cov_fchk_real_regression ${CMAKE_BINARY_DIR}/real_gaussian16_dvb_sp.fchk)
    set_tests_properties(fchk_real_regression PROPERTIES SKIP_RETURN_CODE 77)

    add_executable(cov_fchk_batch_regression tests/fchk_batch_regression.cpp)
    target_link_libraries(cov_fchk_batch_regression PRIVATE
        cov_parser cov_orbitals cov_molecule)

    add_executable(cov_scientific_audit_dump tests/scientific_audit_dump.cpp)
    target_link_libraries(cov_scientific_audit_dump PRIVATE cov_parser cov_molecule cov_orbitals)
    add_executable(cov_pi_diagram_probe tests/pi_diagram_probe.cpp)
    target_link_libraries(cov_pi_diagram_probe PRIVATE cov_parser cov_molecule cov_orbitals)
    add_executable(cov_pi_channel_smoke tests/pi_channel_smoke.cpp)
    target_link_libraries(cov_pi_channel_smoke PRIVATE cov_orbitals)
    add_test(NAME pi_channel_smoke COMMAND cov_pi_channel_smoke)
    add_executable(cov_pi_mode_scope_smoke tests/pi_mode_scope_smoke.cpp)
    target_link_libraries(cov_pi_mode_scope_smoke PRIVATE cov_parser cov_orbitals)
    add_test(NAME pi_mode_scope_smoke COMMAND cov_pi_mode_scope_smoke)
    add_executable(cov_mo_composition_ledger_smoke tests/mo_composition_ledger_smoke.cpp)
    target_link_libraries(cov_mo_composition_ledger_smoke PRIVATE cov_parser cov_orbitals)
    add_test(NAME mo_composition_ledger_smoke COMMAND cov_mo_composition_ledger_smoke)
    add_executable(cov_mo_brief_scope_probe tests/mo_brief_scope_probe.cpp)
    target_link_libraries(cov_mo_brief_scope_probe PRIVATE cov_parser cov_orbitals)

    # These collectors emit raw records for independently frozen control
    # reviewers. A successful process alone is not a scientific verdict.
    foreach(probe symmetry_scope energy_gap_scope pi_prior viewport_layout
                  local_irrep_scope local_projection angular_generators topology_evidence)
        add_executable(cov_${probe}_probe tests/${probe}_probe.cpp)
        target_link_libraries(cov_${probe}_probe PRIVATE cov_orbitals)
        target_compile_options(cov_${probe}_probe PRIVATE
            $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
    endforeach()

    add_executable(cov_reference_manifest_regression
        tests/reference_manifest_regression.cpp)
    target_link_libraries(cov_reference_manifest_regression PRIVATE
        cov_molecule cov_orbitals cov_parser)

    add_executable(cov_gaussian_log_smoke tests/gaussian_log_smoke.cpp)
    target_link_libraries(cov_gaussian_log_smoke PRIVATE cov_parser)
    add_test(NAME gaussian_log_smoke COMMAND cov_gaussian_log_smoke)

    add_executable(cov_wavefunction_analysis_smoke tests/wavefunction_analysis_smoke.cpp)
    target_link_libraries(cov_wavefunction_analysis_smoke PRIVATE cov_parser)
    add_test(NAME wavefunction_analysis_smoke COMMAND cov_wavefunction_analysis_smoke)

    add_executable(cov_density_source_regression tests/density_source_regression.cpp)
    target_link_libraries(cov_density_source_regression PRIVATE cov_parser)
    add_test(NAME density_source_regression COMMAND cov_density_source_regression)

    # Read-only evidence from the same production parser/analysis library.
    add_executable(cov_density_source_probe tests/density_source_probe.cpp)
    target_link_libraries(cov_density_source_probe PRIVATE cov_parser)

    add_executable(cov_overlap_large_smoke tests/overlap_large_smoke.cpp)
    target_link_libraries(cov_overlap_large_smoke PRIVATE cov_parser)
    add_test(NAME overlap_large_smoke COMMAND cov_overlap_large_smoke)
    set_tests_properties(overlap_large_smoke PROPERTIES TIMEOUT 120)

    add_executable(cov_multicentre_classification_smoke tests/multicentre_classification_smoke.cpp)
    target_link_libraries(cov_multicentre_classification_smoke PRIVATE cov_parser)
    add_test(NAME multicentre_classification_smoke COMMAND cov_multicentre_classification_smoke)

    add_executable(cov_chemical_valence_smoke tests/chemical_valence_smoke.cpp)
    target_link_libraries(cov_chemical_valence_smoke PRIVATE cov_parser cov_orbitals)
    add_test(NAME chemical_valence_smoke COMMAND cov_chemical_valence_smoke)

    add_executable(cov_delocalised_pi_smoke tests/delocalised_pi_smoke.cpp)
    target_link_libraries(cov_delocalised_pi_smoke PRIVATE cov_parser cov_orbitals)
    target_compile_options(cov_delocalised_pi_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
    add_test(NAME delocalised_pi_smoke COMMAND cov_delocalised_pi_smoke)

    add_executable(cov_symmetry_smoke tests/symmetry_smoke.cpp)
    target_link_libraries(cov_symmetry_smoke PRIVATE cov_parser)
    add_test(NAME symmetry_smoke COMMAND cov_symmetry_smoke)

    add_executable(cov_point_group_catalog_smoke tests/point_group_catalog_smoke.cpp)
    target_link_libraries(cov_point_group_catalog_smoke PRIVATE cov_parser)
    add_test(NAME point_group_catalog_smoke COMMAND cov_point_group_catalog_smoke)

    add_executable(cov_point_group_irreps_smoke tests/point_group_irreps_smoke.cpp)
    target_link_libraries(cov_point_group_irreps_smoke PRIVATE cov_parser)
    add_test(NAME point_group_irreps_smoke COMMAND cov_point_group_irreps_smoke)

    add_executable(cov_coordination_geometry_smoke tests/coordination_geometry_smoke.cpp)
    target_link_libraries(cov_coordination_geometry_smoke PRIVATE cov_orbitals)
    add_test(NAME coordination_geometry_smoke COMMAND cov_coordination_geometry_smoke)

    add_executable(cov_local_orbital_symmetry_smoke tests/local_orbital_symmetry_smoke.cpp)
    target_link_libraries(cov_local_orbital_symmetry_smoke PRIVATE cov_orbitals)
    add_test(NAME local_orbital_symmetry_smoke COMMAND cov_local_orbital_symmetry_smoke)

    add_executable(cov_local_geometry_smoke tests/local_geometry_smoke.cpp)
    target_link_libraries(cov_local_geometry_smoke PRIVATE cov_orbitals)
    add_test(NAME local_geometry_smoke COMMAND cov_local_geometry_smoke)

    add_executable(cov_orbital_tracking_smoke tests/orbital_tracking_smoke.cpp)
    target_link_libraries(cov_orbital_tracking_smoke PRIVATE cov_orbitals)
    add_test(NAME orbital_tracking_smoke COMMAND cov_orbital_tracking_smoke)

    # Manual external-corpus regression. It is intentionally not registered
    # with CTest because the Gaussian scan fixtures live outside the source
    # tree and are not available on ordinary developer/CI hosts.
    add_executable(cov_scan_manifest_regression
        tests/scan_manifest_regression.cpp)
    target_link_libraries(cov_scan_manifest_regression PRIVATE
        cov_parser cov_orbitals)

    add_executable(cov_pi_topology_smoke tests/pi_topology_smoke.cpp)
    target_link_libraries(cov_pi_topology_smoke PRIVATE cov_orbitals)
    add_test(NAME pi_topology_smoke COMMAND cov_pi_topology_smoke)
    add_executable(cov_topology_graph_smoke tests/topology_graph_smoke.cpp)
    target_link_libraries(cov_topology_graph_smoke PRIVATE cov_parser)
    add_test(NAME topology_graph_smoke COMMAND cov_topology_graph_smoke)

    add_executable(cov_interaction_graph_smoke tests/interaction_graph_smoke.cpp)
    target_link_libraries(cov_interaction_graph_smoke PRIVATE cov_molecule cov_parser)
    add_test(NAME interaction_graph_smoke COMMAND cov_interaction_graph_smoke)

    add_executable(cov_coordination_frame_contract_smoke
        tests/coordination_frame_contract_smoke.cpp)
    target_link_libraries(cov_coordination_frame_contract_smoke PRIVATE cov_orbitals)
    add_test(NAME coordination_frame_contract_smoke
        COMMAND cov_coordination_frame_contract_smoke)

    add_executable(cov_orbital_symmetry_smoke tests/orbital_symmetry_smoke.cpp)
    target_link_libraries(cov_orbital_symmetry_smoke PRIVATE cov_parser)
    add_test(NAME orbital_symmetry_smoke COMMAND cov_orbital_symmetry_smoke)

    add_executable(cov_linear_orbital_symmetry_smoke tests/linear_orbital_symmetry_smoke.cpp)
    target_link_libraries(cov_linear_orbital_symmetry_smoke PRIVATE cov_parser)
    add_test(NAME linear_orbital_symmetry_smoke COMMAND cov_linear_orbital_symmetry_smoke)

    add_executable(cov_orbital_view_smoke tests/orbital_view_smoke.cpp)
    target_link_libraries(cov_orbital_view_smoke PRIVATE cov_orbitals)
    add_test(NAME orbital_view_smoke COMMAND cov_orbital_view_smoke)

    add_executable(cov_mo_diagram_smoke tests/mo_diagram_smoke.cpp)
    target_link_libraries(cov_mo_diagram_smoke PRIVATE cov_orbitals cov_parser cov_molecule)
    add_test(NAME mo_diagram_smoke COMMAND cov_mo_diagram_smoke)

    add_executable(cov_energy_axis_smoke tests/energy_axis_smoke.cpp)
    target_link_libraries(cov_energy_axis_smoke PRIVATE cov_orbitals)
    add_test(NAME energy_axis_smoke COMMAND cov_energy_axis_smoke)

    add_executable(cov_cp_reference_layout_smoke tests/cp_reference_layout_smoke.cpp)
    target_link_libraries(cov_cp_reference_layout_smoke PRIVATE cov_orbitals)
    target_compile_options(cov_cp_reference_layout_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
    add_test(NAME cp_reference_layout_smoke COMMAND cov_cp_reference_layout_smoke)

    add_executable(cov_mo_diagram_export_semantics_smoke
        tests/mo_diagram_export_semantics_smoke.cpp)
    target_link_libraries(cov_mo_diagram_export_semantics_smoke PRIVATE cov_orbitals)
    add_test(NAME mo_diagram_export_semantics_smoke
        COMMAND cov_mo_diagram_export_semantics_smoke)

    add_executable(cov_mo_diagram_snapshot_smoke tests/mo_diagram_snapshot_smoke.cpp)
    add_executable(cov_nbo_layout_probe tests/nbo_layout_probe.cpp)
    target_link_libraries(cov_nbo_layout_probe PRIVATE cov_ui)
    target_compile_options(cov_nbo_layout_probe PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
    target_link_libraries(cov_mo_diagram_snapshot_smoke PRIVATE cov_ui)
    target_compile_options(cov_mo_diagram_snapshot_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
    add_test(NAME mo_diagram_snapshot_smoke COMMAND cov_mo_diagram_snapshot_smoke)
    add_executable(cov_nbo_aomo_labels_smoke tests/nbo_aomo_labels_smoke.cpp)
    add_executable(cov_nbo_mode_scope_probe tests/nbo_mode_scope_probe.cpp)
    target_link_libraries(cov_nbo_mode_scope_probe PRIVATE cov_ui)
    add_executable(cov_nbo_salc_names_probe tests/nbo_salc_names_probe.cpp)
    target_link_libraries(cov_nbo_salc_names_probe PRIVATE cov_ui)
    target_compile_options(cov_nbo_salc_names_probe PRIVATE $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
    add_executable(cov_orbital_runtime_probe tests/orbital_runtime_probe.cpp)
    target_link_libraries(cov_orbital_runtime_probe PRIVATE cov_ui)
    target_compile_options(cov_orbital_runtime_probe PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/utf-8>)
    add_executable(cov_atomic_orbital_symmetry_smoke tests/atomic_orbital_symmetry_smoke.cpp)
    target_link_libraries(cov_atomic_orbital_symmetry_smoke PRIVATE cov_parser)
    target_compile_options(cov_atomic_orbital_symmetry_smoke PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/utf-8>)
    add_test(NAME atomic_orbital_symmetry_smoke COMMAND cov_atomic_orbital_symmetry_smoke)
    add_executable(cov_orbital_symmetry_components_smoke tests/orbital_symmetry_components_smoke.cpp)
    add_executable(cov_orbital_symmetry_components_probe tests/orbital_symmetry_components_probe.cpp)
    target_link_libraries(cov_orbital_symmetry_components_probe PRIVATE cov_ui)
    target_compile_options(cov_orbital_symmetry_components_probe PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/utf-8>)
    target_link_libraries(cov_orbital_symmetry_components_smoke PRIVATE cov_ui)
    target_compile_options(cov_orbital_symmetry_components_smoke PRIVATE $<$<CXX_COMPILER_ID:MSVC>:/utf-8>)
    add_test(NAME orbital_symmetry_components_smoke COMMAND cov_orbital_symmetry_components_smoke)
    target_link_libraries(cov_nbo_aomo_labels_smoke PRIVATE cov_ui)
    target_compile_options(cov_nbo_aomo_labels_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
    add_test(NAME nbo_aomo_labels_smoke COMMAND cov_nbo_aomo_labels_smoke)
    add_executable(cov_nbo_aomo_hover_smoke tests/nbo_aomo_hover_smoke.cpp)
    target_link_libraries(cov_nbo_aomo_hover_smoke PRIVATE cov_ui)
    target_compile_options(cov_nbo_aomo_hover_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
    add_test(NAME nbo_aomo_hover_smoke COMMAND cov_nbo_aomo_hover_smoke)
    add_executable(cov_mo_diagram_layout_smoke tests/mo_diagram_layout_smoke.cpp)
    target_link_libraries(cov_mo_diagram_layout_smoke PRIVATE cov_ui)
    add_test(NAME mo_diagram_layout_smoke COMMAND cov_mo_diagram_layout_smoke)
    add_executable(cov_validation_navigation_smoke tests/validation_navigation_smoke.cpp)
    target_link_libraries(cov_validation_navigation_smoke PRIVATE cov_validation_navigation cov_ui)
    add_test(NAME validation_navigation_smoke COMMAND cov_validation_navigation_smoke)

    add_executable(cov_molecule_style_smoke tests/molecule_style_smoke.cpp)
    target_link_libraries(cov_molecule_style_smoke PRIVATE cov_molecule cov_parser)
    add_test(NAME molecule_style_smoke COMMAND cov_molecule_style_smoke)

    add_executable(cov_ui_smoke tests/ui_smoke.cpp)
    target_link_libraries(cov_ui_smoke PRIVATE cov_ui)
    target_compile_options(cov_ui_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
    add_test(NAME ui_smoke COMMAND cov_ui_smoke)

    add_executable(cov_ui_fchk_smoke tests/ui_fchk_smoke.cpp)
    target_link_libraries(cov_ui_fchk_smoke PRIVATE cov_ui)
    target_compile_options(cov_ui_fchk_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
    add_test(NAME ui_fchk_smoke COMMAND cov_ui_fchk_smoke)

    add_executable(cov_ui_orbital_localisation_smoke
        tests/ui_orbital_localisation_smoke.cpp)
    target_link_libraries(cov_ui_orbital_localisation_smoke PRIVATE cov_ui)
    target_compile_options(cov_ui_orbital_localisation_smoke PRIVATE
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>
    )
    add_test(NAME ui_orbital_localisation_smoke
        COMMAND cov_ui_orbital_localisation_smoke)

    add_executable(cov_ui_details_resize_smoke tests/ui_details_resize_smoke.cpp)
    target_link_libraries(cov_ui_details_resize_smoke PRIVATE cov_ui)
    add_test(NAME ui_details_resize_smoke COMMAND cov_ui_details_resize_smoke)
    add_executable(cov_point_group_geometry_smoke tests/point_group_geometry_smoke.cpp)
    target_link_libraries(cov_point_group_geometry_smoke PRIVATE cov_parser)
    target_compile_definitions(cov_point_group_geometry_smoke PRIVATE COV_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME point_group_geometry_smoke COMMAND cov_point_group_geometry_smoke)
    add_executable(cov_orbital_naming_probe tests/orbital_naming_probe.cpp)
    target_link_libraries(cov_orbital_naming_probe PRIVATE cov_ui)
endif()
