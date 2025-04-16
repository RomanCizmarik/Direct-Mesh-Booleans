cmake_minimum_required(VERSION 3.24)

# For FaRMA to work
# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64
# processors, do not support any of the two.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)


# Project name
project(DirectMeshBooleans)

set(TBB_TEST OFF CACHE BOOL " " FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL " " FORCE)

add_library(DirectMeshBooleans)

target_compile_definitions(DirectMeshBooleans PUBLIC TBB_PARALLEL=1)
target_compile_definitions(DirectMeshBooleans PUBLIC CINOLIB_USES_SHEWCHUK_PREDICATES)

#file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*" "${CMAKE_CURRENT_SOURCE_DIR}/include/*")
file(GLOB_RECURSE ${PROJECT_NAME}_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*" "${CMAKE_CURRENT_SOURCE_DIR}/include/*")
target_sources(${PROJECT_NAME} PRIVATE ${${PROJECT_NAME}_SOURCES})

set(${PROJECT_NAME}_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" PATH)

target_include_directories(DirectMeshBooleans PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

include(FetchContent)

message("Fetch Fast And Robust Mesh Arrangements")
FetchContent_Declare(mesh_arrangement
                    GIT_REPOSITORY "https://github.com/gcherchi/FastAndRobustMeshArrangements.git"
                    GIT_TAG bf7eb71da991a61ff5414946a4b2754bbd327e41) #  Apr 18, 2024 "Merge pull request #23 from gcherchi/dev-mf"
FetchContent_MakeAvailable(mesh_arrangement)
target_include_directories(${PROJECT_NAME} PUBLIC ${mesh_arrangement_SOURCE_DIR}/code/)

message("Fetch OpenMesh")
FetchContent_Declare(OpenMesh
                    GIT_REPOSITORY "https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git"
                    GIT_TAG 24b8690bcd525c66068754d3eba32fc6618f5085) #  Apr 08, 2025 "Merge branch 'fix-retain-deprecation-warning' into 'master'"

# Prevent OpenMesh from installing unnecessary things (optional)
set(OPENMESH_BUILD_UNIT_TESTS OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_DECIMATER OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_VDPROG OFF CACHE BOOL "" FORCE)
set(OPENMESH_BUILD_APPS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(OpenMesh)

#target_link_libraries(DirectMeshBooleans
#  PRIVATE
#  mesh_arrangement
#)

# Compiler-specific options for FaRMA to work
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC "/fp:strict")
	# use intrinsic functions
	target_compile_options(${PROJECT_NAME} PUBLIC "/Oi")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	endif()
	# reserve enough stack size
	target_link_options(${PROJECT_NAME} PUBLIC "/STACK:8421376")
	# turn off annoying warnings
	target_compile_options(${PROJECT_NAME} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
else()
	# set standard optimization level
	target_compile_options(${PROJECT_NAME} PUBLIC -O2)
	# reserve enough stack size	
	target_compile_options(${PROJECT_NAME} PUBLIC -Wl,-z,stacksize=8421376)
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME} PUBLIC -frounding-math)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
	endif()
endif()

#if (MSVC)
#	target_link_options(DirectMeshBooleans PRIVATE "/SUBSYSTEM:CONSOLE")
#endif()

#if (MSVC)
#  add_custom_command(TARGET DirectMeshBooleans POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:DirectMeshBooleans> "$<TARGET_FILE_DIR:DirectMeshBooleans>/"
#    COMMAND mt.exe -manifest \"${CMAKE_CURRENT_SOURCE_DIR}/resources/manifest.xml\" "-inputresource:\"$<TARGET_FILE:DirectMeshBooleans>\"\;#1" "-outputresource:\"$<TARGET_FILE:DirectMeshBooleans>\"\;#1"
#    COMMAND_EXPAND_LISTS
#  )

#  add_custom_command(TARGET DirectMeshBooleans POST_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:DirectMeshBooleans> "$<TARGET_FILE_DIR:DirectMeshBooleans>/"
#    COMMAND_EXPAND_LISTS
#  )
#endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/test")