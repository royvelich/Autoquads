cmake_minimum_required(VERSION 3.5)
# set the project name
project(optimization_lib)

option(OPTIMIZATION_USE_WORHP            "Use worhp"                    OFF)

file(GLOB sources
	src/solvers/*.cpp
	src/objective_functions/*.cpp
	include/solvers/*.h
	include/objective_functions/*.h
	include/utils.h
	)


add_library(${PROJECT_NAME} ${sources})
add_library(rds::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME}
	PRIVATE
		${CMAKE_SOURCE_DIR}
		${PROJECT_SOURCE_DIR}/include
		${CMAKE_SOURCE_DIR}/ext/worhp/worhp)
if(OPTIMIZATION_USE_WORHP)
	add_library(worhp SHARED IMPORTED)
	set_property(TARGET worhp PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/ext/worhp/bin/worhp.dll)
	set_property(TARGET worhp PROPERTY IMPORTED_IMPLIB ${CMAKE_SOURCE_DIR}/ext/worhp/bin/worhp.lib)
endif()

target_link_libraries(${PROJECT_NAME}
	PRIVATE
		igl::core	
		igl::matlab
		)
if(OPTIMIZATION_USE_WORHP)
target_link_libraries(${PROJECT_NAME}
	PRIVATE
	worhp
		)

endif()