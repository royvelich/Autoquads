cmake_minimum_required(VERSION 3.15)
project(console_app)

# Sources
file(GLOB EXTERNAL_SOURCES
	${CMAKE_SOURCE_DIR}/natvis/eigen.natvis)

file(GLOB INTERNAL_SOURCES
	src/main.cpp)

set(SOURCES ${INTERNAL_SOURCES} ${EXTERNAL_SOURCES})

# Add Library Target
add_executable(${PROJECT_NAME} ${SOURCES})

# Include Directories
target_include_directories(${PROJECT_NAME} 
	PRIVATE
		${Boost_INCLUDE_DIRS})

# Link Libraries
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        igl::core
        igl::opengl
        rds::optimization_lib)

find_package(OpenMP)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

# Properties
set_target_properties(${PROJECT_NAME} PROPERTIES
	CXX_STANDARD 17
	VS_GLOBAL_UseIntelMKL "Sequential"
	VS_GLOBAL_UseIntelTBB "Yes")

if (MSVC)
	# Turn on the __cplusplus flag in MSVC, so the __cplusplus macro will report the correct C++ version
	# https://docs.microsoft.com/en-us/cpp/build/reference/zc-cplusplus?view=vs-2019
	# http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1309
	target_compile_options(${PROJECT_NAME} PRIVATE /Zc:__cplusplus)
endif()

# Source Tree
source_group(TREE ${PROJECT_SOURCE_DIR} FILES ${INTERNAL_SOURCES})