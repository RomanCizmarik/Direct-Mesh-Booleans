add_executable(${PROJECT_NAME}Test)

file(GLOB_RECURSE ${PROJECT_NAME}Test_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*" "${CMAKE_CURRENT_SOURCE_DIR}/include/*")
target_sources(${PROJECT_NAME}Test PRIVATE ${${PROJECT_NAME}Test_SOURCES})

target_include_directories(${PROJECT_NAME}Test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(${PROJECT_NAME}Test PRIVATE ${${PROJECT_NAME}_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME}Test PUBLIC OpenMeshCore OpenMeshTools)

if (MSVC)
	target_link_options(${PROJECT_NAME}Test PRIVATE "/SUBSYSTEM:CONSOLE")
endif()