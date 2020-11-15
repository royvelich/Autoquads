cmake_minimum_required(VERSION 3.15)
project(optimization_lib)

# Sources
file(GLOB SOURCES
	src/core/updatable_object.cpp
	src/data_providers/mesh_wrapper.cpp
	src/data_providers/mesh_data_provider.cpp
	src/data_providers/data_provider.cpp
	src/data_providers/empty_data_provider.cpp
	src/data_providers/coordinate_data_provider.cpp
	src/data_providers/coordinate_diff_data_provider.cpp
	src/data_providers/cross_coordinate_diff_data_provider.cpp	
	src/data_providers/plain_data_provider.cpp
	src/data_providers/edge_pair_data_provider.cpp
	src/data_providers/face_fan_data_provider.cpp
	src/data_providers/face_data_provider.cpp
	src/objective_functions/objective_function_base.cpp
	src/objective_functions/objective_function.cpp
	src/objective_functions/concrete_objective.cpp
	src/objective_functions/dense_objective_function.cpp
	src/objective_functions/sparse_objective_function.cpp
	src/objective_functions/summation_objective.cpp
	src/objective_functions/composite_objective.cpp
	src/objective_functions/symmetric_dirichlet_objective.cpp
	src/objective_functions/separation_objective.cpp
	src/objective_functions/periodic_objective.cpp
	src/objective_functions/parabolic_objective.cpp
	src/objective_functions/coordinate_objective.cpp
	src/objective_functions/coordinate_diff_objective.cpp
	src/objective_functions/cross_coordinate_diff_objective.cpp
	src/objective_functions/seamless_objective.cpp
	src/objective_functions/position/patch_position_objective.cpp
	src/objective_functions/position/face_position_objective.cpp
	src/objective_functions/position/face_vertices_position_objective.cpp
	src/objective_functions/position/face_barycenter_position_objective.cpp
	src/objective_functions/edge_pair/edge_pair_objective.cpp
	src/objective_functions/edge_pair/edge_pair_angle_objective.cpp
	src/objective_functions/edge_pair/edge_pair_length_objective.cpp
	src/objective_functions/edge_pair/edge_pair_translation_objective.cpp
	src/objective_functions/edge_pair/edge_pair_integer_translation_objective.cpp
	src/objective_functions/singularity/singular_point_position_objective.cpp
	src/objective_functions/singularity/singular_points_position_objective.cpp
	src/iterative_methods/iterative_method.cpp
	src/iterative_methods/newton_method.cpp
	src/iterative_methods/gradient_descent.cpp
	src/solvers/solver.cpp
	src/solvers/eigen_sparse_solver.cpp
	src/solvers/pardiso_solver.cpp
	include/core/core.h
	include/core/utils.h
	include/core/updatable_object.h
	include/data_providers/mesh_wrapper.h
	include/data_providers/mesh_data_provider.h
	include/data_providers/data_provider.h
	include/data_providers/empty_data_provider.h
	include/data_providers/coordinate_data_provider.h
	include/data_providers/coordinate_diff_data_provider.h
	include/data_providers/cross_coordinate_diff_data_provider.h
	include/data_providers/plain_data_provider.h
	include/data_providers/edge_pair_data_provider.h
	include/data_providers/face_fan_data_provider.h
	include/data_providers/face_data_provider.h
	include/objective_functions/objective_function_base.h	
	include/objective_functions/objective_function.h
	include/objective_functions/concrete_objective.h	
	include/objective_functions/sparse_objective_function.h
	include/objective_functions/dense_objective_function.h
	include/objective_functions/summation_objective.h
	include/objective_functions/composite_objective.h
	include/objective_functions/symmetric_dirichlet_objective.h
	include/objective_functions/separation_objective.h
	include/objective_functions/periodic_objective.h
	include/objective_functions/parabolic_objective.h
	include/objective_functions/coordinate_objective.h
	include/objective_functions/coordinate_diff_objective.h
	include/objective_functions/cross_coordinate_diff_objective.h
	include/objective_functions/seamless_objective.h
	include/objective_functions/position/patch_position_objective.h
	include/objective_functions/position/face_position_objective.h
	include/objective_functions/position/face_vertices_position_objective.h
	include/objective_functions/position/face_barycenter_position_objective.h
	include/objective_functions/edge_pair/edge_pair_objective.h
	include/objective_functions/edge_pair/edge_pair_angle_objective.h
	include/objective_functions/edge_pair/edge_pair_length_objective.h
	include/objective_functions/edge_pair/edge_pair_translation_objective.h
	include/objective_functions/edge_pair/edge_pair_integer_translation_objective.h
	include/objective_functions/singularity/singular_point_position_objective.h
	include/objective_functions/singularity/singular_points_position_objective.h
	include/iterative_methods/iterative_method.h
	include/iterative_methods/newton_method.h
	include/iterative_methods/gradient_descent.h
	include/solvers/solver.h	
	include/solvers/eigen_sparse_solver.h
	include/solvers/pardiso_solver.h)

# Add Library Target
add_library(${PROJECT_NAME} ${SOURCES})
add_library(rds::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

# Include Directories
target_include_directories(${PROJECT_NAME}
	PRIVATE
		${Boost_INCLUDE_DIRS}
		${PROJECT_SOURCE_DIR}/include)

# Link Libraries
target_link_libraries(${PROJECT_NAME}
	PRIVATE
		igl::core)

find_package(OpenMP)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

# Properties
set_target_properties(${PROJECT_NAME} PROPERTIES
	CXX_STANDARD 17
	VS_GLOBAL_UseIntelMKL "Parallel"
	VS_GLOBAL_UseIntelTBB "Yes")

if (MSVC)
	# Turn on the __cplusplus flag in MSVC, so the __cplusplus macro will report the correct C++ version
	# https://docs.microsoft.com/en-us/cpp/build/reference/zc-cplusplus?view=vs-2019
	# http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1309
	target_compile_options(${PROJECT_NAME} PRIVATE /Zc:__cplusplus)
endif()

# Source Tree
source_group(TREE ${PROJECT_SOURCE_DIR} FILES ${SOURCES})