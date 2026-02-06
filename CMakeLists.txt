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
project(DirectMeshBooleans VERSION 1.0.0 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(TBB_TEST OFF CACHE BOOL " " FORCE)
set(TBB_EXAMPLES OFF CACHE BOOL " " FORCE)

add_library(DirectMeshBooleans)

add_compile_definitions(DirectMeshBooleans PUBLIC _HAS_STD_BYTE=0)
target_compile_definitions(DirectMeshBooleans PUBLIC TBB_PARALLEL=1)
target_compile_definitions(DirectMeshBooleans PUBLIC CINOLIB_USES_SHEWCHUK_PREDICATES)
target_compile_definitions(DirectMeshBooleans PUBLIC _HAS_STD_BYTE=0)

set(${PROJECT_NAME}_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" PATH)

#target_include_directories(DirectMeshBooleans PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_include_directories(DirectMeshBooleans PUBLIC  
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  
    $<INSTALL_INTERFACE:include>  # <prefix>/include
)

include(FetchContent)

message("Fetch Eigen")
FetchContent_Declare(Eigen
                    GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
                    GIT_TAG d6689a15d75e0af62a5f0f6112f6b75c50947eea) #  Nov 16, 2024 "Replace instances of EIGEN_CONSTEXPR macro"

set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
set(EIGEN_INSTALL ON)
#set(INCLUDE_INSTALL_DIR ON)
set(EIGEN_BUILD_CMAKE_PACKAGE ON)
FetchContent_MakeAvailable(Eigen)
#target_include_directories(${PROJECT_NAME} PUBLIC ${Eigen_SOURCE_DIR})

#message(${Eigen_SOURCE_DIR})


message("Fetch Fast And Robust Mesh Arrangements")

set(FaRMA_patch git apply ${CMAKE_CURRENT_SOURCE_DIR}/patches/FaRMA_patch.patch)

#FetchContent_Declare(hiberlite
#        GIT_REPOSITORY https://github.com/paulftw/hiberlite
#        GIT_TAG master
#        PATCH_COMMAND ${hiberlite_patch}
#        UPDATE_DISCONNECTED 1
#        )

FetchContent_Declare(mesh_arrangement
                    GIT_REPOSITORY "https://github.com/gcherchi/FastAndRobustMeshArrangements.git"
                    GIT_TAG bf7eb71da991a61ff5414946a4b2754bbd327e41 #  Apr 18, 2024 "Merge pull request #23 from gcherchi/dev-mf"
					PATCH_COMMAND ${FaRMA_patch}
					UPDATE_DISCONNECTED 1)


FetchContent_MakeAvailable(mesh_arrangement)
FetchContent_MakeAvailable(indirect_predicates)
FetchContent_MakeAvailable(cinolib)

#target_include_directories(DirectMeshBooleans PUBLIC  
#    $<BUILD_INTERFACE:${mesh_arrangement_SOURCE_DIR}/code/>  
#    $<INSTALL_INTERFACE:BUILD/_deps/mesh_arrangement-src/code/>  # <prefix>/include
#)

target_include_directories(${PROJECT_NAME} PUBLIC ${mesh_arrangement_SOURCE_DIR}/code/)
target_include_directories(${PROJECT_NAME} PUBLIC ${mesh_arrangement_SOURCE_DIR}/external/oneTBB/include)
target_include_directories(${PROJECT_NAME} PUBLIC ${mesh_arrangement_SOURCE_DIR}/external/abseil-cpp/)
target_include_directories(${PROJECT_NAME} PUBLIC ${mesh_arrangement_SOURCE_DIR}/external/parallel-hashmap/)
target_include_directories(${PROJECT_NAME} PUBLIC ${mesh_arrangement_SOURCE_DIR}/)
target_include_directories(${PROJECT_NAME} PUBLIC ${cinolib_SOURCE_DIR}/include/)
target_include_directories(${PROJECT_NAME} PUBLIC ${indirect_predicates_SOURCE_DIR}/include/)



#target_include_directories(${PROJECT_NAME} PUBLIC ${FETCHCONTENT_BASE_DIR}/mesh_arrangement-src/code/)
#target_include_directories(${PROJECT_NAME} PUBLIC ${FETCHCONTENT_BASE_DIR}/mesh_arrangement-src/external/oneTBB/include)
#target_include_directories(${PROJECT_NAME} PUBLIC ${FETCHCONTENT_BASE_DIR}/mesh_arrangement-src/external/abseil-cpp/)
#target_include_directories(${PROJECT_NAME} PUBLIC ${FETCHCONTENT_BASE_DIR}/cinolib-src/include)
#target_include_directories(${PROJECT_NAME} PUBLIC ${FETCHCONTENT_BASE_DIR}/indirect_predicates-src/include)

#message(${indirect_predicates_SOURCE_DIR}/include/)
#message(${mesh_arrangement_SOURCE_DIR}/code/)
#message(${FETCHCONTENT_BASE_DIR}/mesh_arrangement-src/code)
#message(${FETCHCONTENT_BASE_DIR}/cinolib-src/include)
#message(${FETCHCONTENT_BASE_DIR}/indirect_predicates-src/include)

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

#file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*" "${CMAKE_CURRENT_SOURCE_DIR}/include/*")
file(GLOB_RECURSE ${PROJECT_NAME}_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*" "${CMAKE_CURRENT_SOURCE_DIR}/include/*" ) #"${mesh_arrangement_SOURCE_DIR}/code/*"
target_sources(${PROJECT_NAME} PRIVATE ${${PROJECT_NAME}_SOURCES})

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

target_link_libraries(${PROJECT_NAME} PUBLIC OpenMeshCore OpenMeshTools eigen shewchuk_predicates tbb) #Eigen3::Eigen
#shewchuk_predicates tbb


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


# Install headers
#install(DIRECTORY include/ DESTINATION include)

# Export target
#install(TARGETS DirectMeshBooleans
#    EXPORT DirectMeshBooleansTargets
#)

# Generate and install CMake package config files
#include(CMakePackageConfigHelpers)

#write_basic_package_version_file(
#    "${CMAKE_CURRENT_BINARY_DIR}/DirectMeshBooleansConfigVersion.cmake"
#    VERSION ${PROJECT_VERSION}
#    COMPATIBILITY SameMajorVersion  
#)

#install(EXPORT DirectMeshBooleansTargets
#    FILE DirectMeshBooleansTargets.cmake
#    NAMESPACE DirectMeshBooleans::
#    DESTINATION lib/cmake/DirectMeshBooleans
#)

#configure_package_config_file(
#    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DirectMeshBooleansConfig.cmake.in"
#    "${CMAKE_CURRENT_BINARY_DIR}/DirectMeshBooleansConfig.cmake"
#    INSTALL_DESTINATION lib/cmake/DirectMeshBooleans
#)

#install(FILES
#    "${CMAKE_CURRENT_BINARY_DIR}/DirectMeshBooleansConfig.cmake"
#    "${CMAKE_CURRENT_BINARY_DIR}/DirectMeshBooleansConfigVersion.cmake"
#    DESTINATION lib/cmake/DirectMeshBooleans
#)


#Exported alias name
#add_library(DirectMeshBooleans::DirectMeshBooleans ALIAS DirectMeshBooleans)


# Build unit tests if requested
option(DMB_BUILD_UNIT_TESTS "Build unit tests" OFF)
if(DMB_BUILD_UNIT_TESTS)
    add_subdirectory(tests)
endif()


# Add example dir if requested
option(DMB_BUILD_EXAMPLES "Build example applications" OFF)
if(DMB_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()