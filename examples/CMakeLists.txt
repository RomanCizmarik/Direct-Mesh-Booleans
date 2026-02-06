#TODO: rename to DirectMeshBooleansExample_Example

#project(DirectMeshBooleansExample LANGUAGES CXX)

# Only use find_package if the project is not being built as part of the same CMake project
#if (NOT TARGET DirectMeshBooleans::DirectMeshBooleans)
#    find_package(DirectMeshBooleans REQUIRED)
#endif()

add_executable(DirectMeshBooleansExample)

file(GLOB_RECURSE DirectMeshBooleansExample_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*" "${CMAKE_CURRENT_SOURCE_DIR}/include/*")
target_sources(DirectMeshBooleansExample PRIVATE ${DirectMeshBooleansExample_SOURCES})


target_include_directories(DirectMeshBooleansExample PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
#target_include_directories(DirectMeshBooleansExampleTest PRIVATE ${DirectMeshBooleansExample_INCLUDE_DIR})
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${mesh_arrangement_SOURCE_DIR}/code/)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${mesh_arrangement_SOURCE_DIR}/external/oneTBB/include)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${mesh_arrangement_SOURCE_DIR}/external/abseil-cpp/)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${cinolib_SOURCE_DIR}/include/)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${indirect_predicates_SOURCE_DIR}/include/)

target_link_libraries(DirectMeshBooleansExample PUBLIC DirectMeshBooleans)  # OpenMeshCore OpenMeshTools Eigen3::Eigen tbb

if (MSVC)
	target_link_options(DirectMeshBooleansExample PRIVATE "/SUBSYSTEM:CONSOLE")
endif()