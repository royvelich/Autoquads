cmake_minimum_required(VERSION 3.15)
project(node_addon)

# Sources
file(GLOB EXTERNAL_SOURCES
	${CMAKE_SOURCE_DIR}/natvis/eigen.natvis)

file(GLOB INTERNAL_SOURCES
	src/addon.cpp
	src/engine.cpp
	include/engine.h)

set(SOURCES ${INTERNAL_SOURCES} ${EXTERNAL_SOURCES})

# Include N-API wrappers
execute_process(COMMAND node -p "require('node-addon-api').include"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE NODE_ADDON_API_DIR)

string(REPLACE "\n" "" NODE_ADDON_API_DIR ${NODE_ADDON_API_DIR})
string(REPLACE "\"" "" NODE_ADDON_API_DIR ${NODE_ADDON_API_DIR})

# Add Library Target
add_library(${PROJECT_NAME} SHARED ${SOURCES} ${CMAKE_JS_SRC})

# Include Directories
target_include_directories(${PROJECT_NAME} 
	PRIVATE
		${CMAKE_JS_INC}
		${CMAKE_SOURCE_DIR}
		${NODE_ADDON_API_DIR}
		${Boost_INCLUDE_DIRS})

# Link Libraries
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        igl::core
        igl::opengl
        rds::optimization_lib
        ${CMAKE_JS_LIB})

find_package(OpenMP)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

# Properties
set_target_properties(${PROJECT_NAME} PROPERTIES
	CXX_STANDARD 17
	PREFIX ""
	SUFFIX ".node"
	VS_GLOBAL_UseIntelMKL "Parallel"
	VS_GLOBAL_UseIntelTBB "Yes")

if (MSVC)
	# Turn on the __cplusplus flag in MSVC, so the __cplusplus macro will report the correct C++ version
	# https://docs.microsoft.com/en-us/cpp/build/reference/zc-cplusplus?view=vs-2019
	# http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1309
	target_compile_options(${PROJECT_NAME} PRIVATE /Zc:__cplusplus)
endif()

# Source Tree
source_group(TREE ${PROJECT_SOURCE_DIR} FILES ${INTERNAL_SOURCES})

# Custom Commands
add_custom_command(TARGET ${PROJECT_NAME}
                   POST_BUILD
                   COMMAND call $(ProjectDir)..\\..\\scripts\\post_build.bat ARGS "$(TargetDir)$(TargetFileName)" "$(ProjectDir)latest_builds\\"
                   COMMENT "Copying .node file to latest_builds folder")