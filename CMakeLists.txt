cmake_minimum_required(VERSION 3.24)

# ============================================================
# Project
# ============================================================
project(DirectMeshBooleans
    VERSION 1.0.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================
# Options
# ============================================================
option(DMB_PROVIDE_OWN_DEPENDENCIES
    "Use externally provided dependencies via find_package"
    OFF
)

option(DMB_BUILD_UNIT_TESTS "Build unit tests" OFF)
option(DMB_BUILD_EXAMPLES   "Build example applications" OFF)

# SIMD options 
# For FaRMA to work
# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64
# processors, do not support any of the two.
set(ENABLE_SSE2 ON)
set(ENABLE_AVX2 ON)

# ============================================================
# Core library
# ============================================================
add_library(DirectMeshBooleans)

target_include_directories(DirectMeshBooleans
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_definitions(DirectMeshBooleans
    PUBLIC
        _HAS_STD_BYTE=0
)

# ============================================================
# Dependencies
# ============================================================

if(DMB_PROVIDE_OWN_DEPENDENCIES)

    message(STATUS "DirectMeshBooleans: using externally provided dependencies")

    # ---- Eigen (ALWAYS REQUIRED) ----
    find_package(Eigen3 CONFIG REQUIRED)

    # ---- OpenMesh ----
    find_package(OpenMesh CONFIG REQUIRED)
    add_library(OpenMeshCore  ALIAS OpenMesh::Core)
    add_library(OpenMeshTools ALIAS OpenMesh::Tools)

    # ---- FaRMA ----
    # External projects MUST provide a FaRMA target
    find_package(FaRMA CONFIG REQUIRED)

else()

    message(STATUS "DirectMeshBooleans: fetching dependencies")

    include(FetchContent)

    # ------------------------------------------------------------
    # Eigen
    # ------------------------------------------------------------
    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG d6689a15d75e0af62a5f0f6112f6b75c50947eea
    )

    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(EIGEN_MPL2_ONLY ON CACHE BOOL "" FORCE)
    set(EIGEN_INSTALL OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(Eigen)

    # ------------------------------------------------------------
    # FaRMA
    # ------------------------------------------------------------
    set(TBB_TEST OFF CACHE BOOL "" FORCE)
    set(TBB_EXAMPLES OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        mesh_arrangement
        GIT_REPOSITORY https://github.com/gcherchi/FastAndRobustMeshArrangements.git
        GIT_TAG bf7eb71da991a61ff5414946a4b2754bbd327e41
        PATCH_COMMAND
            git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/FaRMA_patch.patch
        UPDATE_DISCONNECTED 1
    )

    FetchContent_MakeAvailable(mesh_arrangement)
    FetchContent_MakeAvailable(indirect_predicates)
    FetchContent_MakeAvailable(cinolib)

    # ------------------------------------------------------------
    # Cinolib (override vendored Eigen!)
    # ------------------------------------------------------------
    add_library(Cinolib INTERFACE)

    target_include_directories(Cinolib
        INTERFACE
            ${cinolib_SOURCE_DIR}/include
    )

    target_link_libraries(Cinolib
        INTERFACE
            Eigen3::Eigen
    )

    target_compile_definitions(Cinolib
        INTERFACE
            CINOLIB_USES_SHEWCHUK_PREDICATES
    )

    # ------------------------------------------------------------
    # FaRMA wrapper target
    # ------------------------------------------------------------
    add_library(FaRMA INTERFACE)

    target_include_directories(FaRMA
        INTERFACE
            ${mesh_arrangement_SOURCE_DIR}/code
            ${mesh_arrangement_SOURCE_DIR}
            ${mesh_arrangement_SOURCE_DIR}/external/abseil-cpp
            ${mesh_arrangement_SOURCE_DIR}/external/parallel-hashmap
            ${indirect_predicates_SOURCE_DIR}/include
    )

    target_compile_definitions(FaRMA
        INTERFACE
            TBB_PARALLEL=1
    )

    target_link_libraries(FaRMA
        INTERFACE
            Cinolib
            shewchuk_predicates
            tbb
    )

    # ------------------------------------------------------------
    # OpenMesh
    # ------------------------------------------------------------
    FetchContent_Declare(
        OpenMesh
        GIT_REPOSITORY https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
        GIT_TAG 24b8690bcd525c66068754d3eba32fc6618f5085
    )

    set(OPENMESH_BUILD_UNIT_TESTS OFF CACHE BOOL "" FORCE)
    set(OPENMESH_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
    set(OPENMESH_BUILD_APPS OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(OpenMesh)

endif()

# ============================================================
# Link DMB (ONLY its real dependencies)
# ============================================================
target_link_libraries(DirectMeshBooleans
    PUBLIC
        FaRMA
        Eigen3::Eigen
        OpenMeshCore
        OpenMeshTools
)

# ============================================================
# Sources
# ============================================================
file(GLOB_RECURSE DMB_SOURCES
    CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*
)

target_sources(DirectMeshBooleans PRIVATE ${DMB_SOURCES})

# ============================================================
# Compiler / platform options (FaRMA requirements)
# ============================================================
if(MSVC)
    target_compile_options(DirectMeshBooleans PUBLIC /fp:strict /Oi)
    if(ENABLE_AVX2)
        target_compile_options(DirectMeshBooleans PUBLIC /arch:AVX2)
    elseif(ENABLE_SSE2)
        target_compile_options(DirectMeshBooleans PUBLIC /arch:SSE2)
    endif()
    target_link_options(DirectMeshBooleans PUBLIC /STACK:8421376)
    target_compile_definitions(DirectMeshBooleans PUBLIC _CRT_SECURE_NO_WARNINGS)
else()
    target_compile_options(DirectMeshBooleans PUBLIC -O2 -frounding-math)
    if(ENABLE_AVX2)
        target_compile_options(DirectMeshBooleans PUBLIC -mavx2)
    elseif(ENABLE_SSE2)
        target_compile_options(DirectMeshBooleans PUBLIC -msse2)
    endif()
endif()

# ============================================================
# Optional components
# ============================================================
if(DMB_BUILD_UNIT_TESTS)
    add_subdirectory(tests)
endif()

if(DMB_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
