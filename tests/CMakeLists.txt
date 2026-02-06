add_executable(DirectMeshBooleansTests)

file(GLOB_RECURSE DirectMeshBooleansTests_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*")
target_sources(DirectMeshBooleansTests PRIVATE ${DirectMeshBooleansTests_SOURCES})

add_compile_definitions(DirectMeshBooleansTests PUBLIC DMB_DATA_PATH=${CMAKE_SOURCE_DIR}/data/)

message("Fetch GTest")
FetchContent_Declare(
                    googletest
                    GIT_REPOSITORY "https://github.com/google/googletest.git"
                    GIT_TAG 52eb8108c5bdec04579160ae17225d66034bd723) #  Apr 30, 2025 "Prepare for v1.17.0"
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)



#target_include_directories(DirectMeshBooleansExample PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
#target_include_directories(DirectMeshBooleansExampleTest PRIVATE ${DirectMeshBooleansExample_INCLUDE_DIR})
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${mesh_arrangement_SOURCE_DIR}/code/)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${mesh_arrangement_SOURCE_DIR}/external/oneTBB/include)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${mesh_arrangement_SOURCE_DIR}/external/abseil-cpp/)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${cinolib_SOURCE_DIR}/include/)
#target_include_directories(DirectMeshBooleansExampleTest PUBLIC ${indirect_predicates_SOURCE_DIR}/include/)

target_link_libraries(DirectMeshBooleansTests PUBLIC DirectMeshBooleans GTest::gtest_main)  # OpenMeshCore OpenMeshTools Eigen3::Eigen tbb

if (MSVC)
	target_link_options(DirectMeshBooleansTests PRIVATE "/SUBSYSTEM:CONSOLE")
endif()