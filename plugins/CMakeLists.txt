cmake_minimum_required(VERSION 3.1)
# set the project name
project(plugins)

# set variable value
set(SOURCES
	deformation_plugin/src/deformation_plugin.cpp
	deformation_plugin/include/deformation_plugin.h
	deformation_plugin/include/app_utils.h)

add_library(${PROJECT_NAME} ${SOURCES})
add_library(rds::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME}
	PRIVATE
		${CMAKE_SOURCE_DIR}
		${PROJECT_SOURCE_DIR}
		${CMAKE_SOURCE_DIR}/libs/
		)

target_link_libraries(${PROJECT_NAME}
	PRIVATE
		igl::core
		igl::opengl
		igl::opengl_glfw
		igl::opengl_glfw_imgui
		igl::matlab
		rds::optimization_lib)		
		  