#pragma once
#ifndef OPTIMIZATION_LIB_EDGE_PAIR_ANGLE_OBJECTIVE_H
#define OPTIMIZATION_LIB_EDGE_PAIR_ANGLE_OBJECTIVE_H

// STL includes
#include <utility>
#include <cmath>
#include <algorithm>
#include <unordered_map>

#define _USE_MATH_DEFINES
#include <cmath>

// Eigen includes
#include <Eigen/Eigenvalues> 

// Optimization lib includes
#include "./periodic_objective.h"
#include "../utils/utils.h"

template<Eigen::StorageOptions StorageOrder_>
class EdgePairAngleObjective : public PeriodicObjective<StorageOrder_>
{
public:
	/**
	 * Public type definitions
	 */
	enum class Properties : uint32_t
	{
		Period = SparseObjectiveFunction<StorageOrder_>::Properties::Count_
	};

	/**
	 * Constructors and destructor
	 */
	EdgePairAngleObjective(const std::shared_ptr<ObjectiveFunctionDataProvider>& objective_function_data_provider, const std::pair<uint64_t, uint64_t>& edge1_indices, const std::pair<uint64_t, uint64_t>& edge2_indices) :
		PeriodicObjective(objective_function_data_provider, "Edge Pair Angle Objective", M_PI / 2),
		edge1_indices_(edge1_indices),
		edge2_indices_(edge2_indices)
	{
		this->Initialize();

		uint64_t edge1_v1_x_index = edge1_indices_.first;
		uint64_t edge1_v1_y_index = edge1_indices_.first + this->image_vertices_count_;
		uint64_t edge1_v2_x_index = edge1_indices_.second;
		uint64_t edge1_v2_y_index = edge1_indices_.second + this->image_vertices_count_;

		uint64_t edge2_v1_x_index = edge2_indices_.first;
		uint64_t edge2_v1_y_index = edge2_indices_.first + this->image_vertices_count_;
		uint64_t edge2_v2_x_index = edge2_indices_.second;
		uint64_t edge2_v2_y_index = edge2_indices_.second + this->image_vertices_count_;

		partial_to_sparse_index_map_.insert({ Partial::E1_V1_X, edge1_v1_x_index });
		partial_to_sparse_index_map_.insert({ Partial::E1_V1_Y, edge1_v1_y_index });
		partial_to_sparse_index_map_.insert({ Partial::E1_V2_X, edge1_v2_x_index });
		partial_to_sparse_index_map_.insert({ Partial::E1_V2_Y, edge1_v2_y_index });
		partial_to_sparse_index_map_.insert({ Partial::E2_V1_X, edge2_v1_x_index });
		partial_to_sparse_index_map_.insert({ Partial::E2_V1_Y, edge2_v1_y_index });
		partial_to_sparse_index_map_.insert({ Partial::E2_V2_X, edge2_v2_x_index });
		partial_to_sparse_index_map_.insert({ Partial::E2_V2_Y, edge2_v2_y_index });

		sparse_index_to_partial_map_.insert({ edge1_v1_x_index, Partial::E1_V1_X });
		sparse_index_to_partial_map_.insert({ edge1_v1_y_index, Partial::E1_V1_Y });
		sparse_index_to_partial_map_.insert({ edge1_v2_x_index, Partial::E1_V2_X });
		sparse_index_to_partial_map_.insert({ edge1_v2_y_index, Partial::E1_V2_Y });
		sparse_index_to_partial_map_.insert({ edge2_v1_x_index, Partial::E2_V1_X });
		sparse_index_to_partial_map_.insert({ edge2_v1_y_index, Partial::E2_V1_Y });
		sparse_index_to_partial_map_.insert({ edge2_v2_x_index, Partial::E2_V2_X });
		sparse_index_to_partial_map_.insert({ edge2_v2_y_index, Partial::E2_V2_Y });
	
		indices_.resize(8);
		indices_.push_back(edge1_v1_x_index);
		indices_.push_back(edge1_v1_y_index);
		indices_.push_back(edge1_v2_x_index);
		indices_.push_back(edge1_v2_y_index);
		indices_.push_back(edge2_v1_x_index);
		indices_.push_back(edge2_v1_y_index);
		indices_.push_back(edge2_v2_x_index);
		indices_.push_back(edge2_v2_y_index);
		std::sort(indices_.begin(), indices_.end());

		for(std::size_t i = 0; i < indices_.size(); i++)
		{
			partial_to_dense_index_map_.insert({ sparse_index_to_partial_map_[indices_.at(i)], i });
			dense_index_to_sparse_index_map_.insert({ i, indices_.at(i) });
		}

		partial_to_first_derivative_sign_map_.insert({ Partial::E1_V1_X,  1 });
		partial_to_first_derivative_sign_map_.insert({ Partial::E1_V1_Y, -1 });
		partial_to_first_derivative_sign_map_.insert({ Partial::E1_V2_X, -1 });
		partial_to_first_derivative_sign_map_.insert({ Partial::E1_V2_Y,  1 });
		partial_to_first_derivative_sign_map_.insert({ Partial::E2_V1_X, -1 });
		partial_to_first_derivative_sign_map_.insert({ Partial::E2_V1_Y,  1 });
		partial_to_first_derivative_sign_map_.insert({ Partial::E2_V2_X,  1 });
		partial_to_first_derivative_sign_map_.insert({ Partial::E2_V2_Y, -1 });

		partial_to_partial_coordinate_type_.insert({ Partial::E1_V1_X, PartialCoordinateType::E1_X });
		partial_to_partial_coordinate_type_.insert({ Partial::E1_V1_Y, PartialCoordinateType::E1_Y });
		partial_to_partial_coordinate_type_.insert({ Partial::E1_V2_X, PartialCoordinateType::E1_X });
		partial_to_partial_coordinate_type_.insert({ Partial::E1_V2_Y, PartialCoordinateType::E1_Y });
		partial_to_partial_coordinate_type_.insert({ Partial::E2_V1_X, PartialCoordinateType::E2_X });
		partial_to_partial_coordinate_type_.insert({ Partial::E2_V1_Y, PartialCoordinateType::E2_Y });
		partial_to_partial_coordinate_type_.insert({ Partial::E2_V2_X, PartialCoordinateType::E2_X });
		partial_to_partial_coordinate_type_.insert({ Partial::E2_V2_Y, PartialCoordinateType::E2_Y });

		for (int i = 0; i < 8; i++)
		{
			auto partial = static_cast<Partial>(i);
			partial_to_first_derivative_value_map_.insert({ partial, 0 });
		}
		
		for(int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				PartialCoordinateType first_partial_type = static_cast<PartialCoordinateType>(i);
				Partial second_partial = static_cast<Partial>(j);
				partial_to_second_derivative_value_map_.insert({ std::make_pair(first_partial_type, second_partial), 0 });
			}
		}
	}

	virtual ~EdgePairAngleObjective()
	{

	}

protected:
	/**
	 * Protected overrides
	 */
	void CalculateValueInner(double& f) override
	{
		f = std::atan2(edge1_y_diff_, edge1_x_diff_) - atan2(edge2_y_diff_, edge2_x_diff_);
	}
	
	void CalculateGradientInner(Eigen::SparseVector<double>& g) override
	{
		for(int i = 0; i < 8; i++)
		{
			Partial partial = static_cast<Partial>(i);
			g.coeffRef(partial_to_sparse_index_map_[partial]) = CalculateFirstPartialDerivative(partial);
		}
	}
	
	void CalculateTripletsInner(std::vector<Eigen::Triplet<double>>& triplets) override
	{
		Eigen::MatrixXd H_dense(8, 8);
		for (int i = 0; i < 8; i++)
		{
			Partial first_partial = static_cast<Partial>(i);
			for (auto j = 0; j < 8; j++)
			{
				Partial second_partial = static_cast<Partial>(j);
				H_dense.coeffRef(partial_to_dense_index_map_[first_partial], partial_to_dense_index_map_[second_partial]) = CalculateSecondPartialDerivative(first_partial, second_partial);;
			}
		}

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H_dense);
		Eigen::MatrixXd D = solver.eigenvalues().asDiagonal();
		Eigen::MatrixXd V = solver.eigenvectors();
		for (auto i = 0; i < 8; i++)
		{
			auto& value = D.coeffRef(i, i);
			if (value < 0)
			{
				value = 0;
			}
		}

		H_dense = V * D * V.inverse();

		auto triplet_index = 0;
		for (auto column = 0; column < 8; column++)
		{
			for (auto row = 0; row <= column; row++)
			{
				const_cast<double&>(triplets[triplet_index].value()) = H_dense(row, column);
				triplet_index++;
			}
		}
	}

	void PreUpdate(const Eigen::VectorXd& x) override
	{
		edge1_v1_ = Eigen::Vector2d(x(edge1_v1_x_index_), x(edge1_v1_y_index_));
		edge1_v2_ = Eigen::Vector2d(x(edge1_v2_x_index_), x(edge1_v2_y_index_));
		
		edge2_v1_ = Eigen::Vector2d(x(edge2_v1_x_index_), x(edge2_v1_y_index_));
		edge2_v2_ = Eigen::Vector2d(x(edge2_v2_x_index_), x(edge2_v2_y_index_));

		edge1_x_diff_ = edge1_v2_.x() - edge1_v1_.x();
		edge1_y_diff_ = edge1_v2_.y() - edge1_v1_.y();
		edge2_x_diff_ = edge2_v2_.x() - edge2_v1_.x();
		edge2_y_diff_ = edge2_v2_.y() - edge2_v1_.y();

		edge1_x_diff_squared_ = edge1_x_diff_ * edge1_x_diff_;
		edge1_y_diff_squared_ = edge1_y_diff_ * edge1_y_diff_;
		edge2_x_diff_squared_ = edge2_x_diff_ * edge2_x_diff_;
		edge2_y_diff_squared_ = edge2_y_diff_ * edge2_y_diff_;

		edge1_squared_norm_ = edge1_x_diff_squared_ + edge1_y_diff_squared_;
		edge2_squared_norm_ = edge2_x_diff_squared_ + edge2_y_diff_squared_;
		
		edge1_quadrupled_norm_ = edge1_squared_norm_ * edge1_squared_norm_;
		edge2_quadrupled_norm_ = edge2_squared_norm_ * edge2_squared_norm_;

		double e1_y_to_e1_squared_norm = edge1_y_diff_ / edge1_squared_norm_;
		double e1_x_to_e1_squared_norm = edge1_x_diff_ / edge1_squared_norm_;
		double e2_y_to_e2_squared_norm = edge2_y_diff_ / edge2_squared_norm_;
		double e2_x_to_e2_squared_norm = edge2_x_diff_ / edge2_squared_norm_;

		partial_to_first_derivative_value_map_[Partial::E1_V1_X] = e1_y_to_e1_squared_norm;
		partial_to_first_derivative_value_map_[Partial::E1_V1_Y] = e1_x_to_e1_squared_norm;
		partial_to_first_derivative_value_map_[Partial::E1_V2_X] = e1_y_to_e1_squared_norm;
		partial_to_first_derivative_value_map_[Partial::E1_V2_Y] = e1_x_to_e1_squared_norm;
		partial_to_first_derivative_value_map_[Partial::E2_V1_X] = e2_y_to_e2_squared_norm;
		partial_to_first_derivative_value_map_[Partial::E2_V1_Y] = e2_x_to_e2_squared_norm;
		partial_to_first_derivative_value_map_[Partial::E2_V2_X] = e2_y_to_e2_squared_norm;
		partial_to_first_derivative_value_map_[Partial::E2_V2_Y] = e2_x_to_e2_squared_norm;

		double e1_diff_prod_to_quad_norm = (2 * edge1_x_diff_ * edge1_y_diff_) / edge1_quadrupled_norm_;
		double e2_diff_prod_to_quad_norm = (2 * edge2_x_diff_ * edge2_y_diff_) / edge2_quadrupled_norm_;
		double e1_squares_diff_prod_to_quad_norm = (edge1_x_diff_squared_ - edge1_y_diff_squared_) / edge1_quadrupled_norm_;
		double e2_squares_diff_prod_to_quad_norm = (edge2_x_diff_squared_ - edge2_y_diff_squared_) / edge2_quadrupled_norm_;

		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E1_V1_X }] = e1_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E1_V1_Y }] = -e1_squares_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E1_V2_X }] = -e1_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E1_V2_Y }] = e1_squares_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E2_V1_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E2_V1_Y }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E2_V2_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_X, Partial::E2_V2_Y }] = 0;

		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E1_V1_X }] = e1_squares_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E1_V1_Y }] = e1_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E1_V2_X }] = -e1_squares_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E1_V2_Y }] = -e1_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E2_V1_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E2_V1_Y }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E2_V2_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E1_Y, Partial::E2_V2_Y }] = 0;

		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E1_V1_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E1_V1_Y }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E1_V2_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E1_V2_Y }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E2_V1_X }] = e2_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E2_V1_Y }] = -e2_squares_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E2_V2_X }] = -e2_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_X, Partial::E2_V2_Y }] = e2_squares_diff_prod_to_quad_norm;

		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E1_V1_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E1_V1_Y }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E1_V2_X }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E1_V2_Y }] = 0;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E2_V1_X }] = e2_squares_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E2_V1_Y }] = e2_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E2_V2_X }] = -e2_squares_diff_prod_to_quad_norm;
		partial_to_second_derivative_value_map_[{ PartialCoordinateType::E2_Y, Partial::E2_V2_Y }] = -e2_diff_prod_to_quad_norm;
	}

private:
	/**
	 * Private type declarations
	 */
	enum class Partial
	{
		E1_V1_X,
		E1_V1_Y,
		E1_V2_X,
		E1_V2_Y,
		E2_V1_X,
		E2_V1_Y,
		E2_V2_X,
		E2_V2_Y
	};

	enum class PartialCoordinateType
	{
		E1_X,
		E1_Y,
		E2_X,
		E2_Y
	};

	/**
	 * Private overrides
	 */
	void CalculateValuePerVertex(Eigen::SparseVector<double>& f_per_vertex) override
	{
		double value = this->GetValue();
		f_per_vertex.coeffRef(edge1_indices_.first) += value;
		f_per_vertex.coeffRef(edge1_indices_.second) += value;
		f_per_vertex.coeffRef(edge2_indices_.first) += value;
		f_per_vertex.coeffRef(edge2_indices_.second) += value;
	}
	
	/**
	 * Private methods
	 */
	double CalculateFirstPartialDerivative(const Partial partial)
	{
		return partial_to_first_derivative_sign_map_[partial] * partial_to_first_derivative_value_map_[partial];
	}

	
	double CalculateSecondPartialDerivative(const Partial first_partial, const Partial second_partial)
	{
		return partial_to_first_derivative_sign_map_[first_partial] * partial_to_second_derivative_value_map_[{ partial_to_partial_coordinate_type_[first_partial], second_partial }];
	}
	
	/**
	 * Private overrides
	 */
	void InitializeTriplets(std::vector<Eigen::Triplet<double>>& triplets) override
	{
		triplets.resize(36);
		auto triplet_index = 0;
		for (auto column = 0; column < 8; column++)
		{
			for (auto row = 0; row <= column; row++)
			{
				triplets[triplet_index] = Eigen::Triplet<double>(indices_.at(row), indices_.at(column), 0);
				triplet_index++;
			}
		}
	}
	
	/**
	 * Private fields
	 */
	std::pair<uint64_t, uint64_t> edge1_indices_;
	std::pair<uint64_t, uint64_t> edge2_indices_;

	std::vector<uint64_t> indices_;
	
	uint64_t edge1_v1_x_index_;
	uint64_t edge1_v1_y_index_;
	uint64_t edge1_v2_x_index_;
	uint64_t edge1_v2_y_index_;

	uint64_t edge2_v1_x_index_;
	uint64_t edge2_v1_y_index_;
	uint64_t edge2_v2_x_index_;
	uint64_t edge2_v2_y_index_;

	Eigen::Vector2d edge1_v1_;
	Eigen::Vector2d edge1_v2_;
	Eigen::Vector2d edge2_v1_;
	Eigen::Vector2d edge2_v2_;

	double edge1_x_diff_;
	double edge1_y_diff_;
	double edge2_x_diff_;
	double edge2_y_diff_;

	double edge1_x_diff_squared_;
	double edge1_y_diff_squared_;
	double edge2_x_diff_squared_;
	double edge2_y_diff_squared_;
	
	double edge1_atan2_dx_;
	double edge1_atan2_dy_;
	double edge2_atan2_dx_;
	double edge2_atan2_dy_;

	double edge1_squared_norm_;
	double edge2_squared_norm_;

	double edge1_quadrupled_norm_;
	double edge2_quadrupled_norm_;

	std::unordered_map<Partial, PartialCoordinateType> partial_to_partial_coordinate_type_;
	std::unordered_map<Partial, uint64_t> partial_to_sparse_index_map_;
	std::unordered_map<uint64_t, Partial> sparse_index_to_partial_map_;
	std::unordered_map<Partial, uint64_t> partial_to_dense_index_map_;
	std::unordered_map<uint64_t, uint64_t> dense_index_to_sparse_index_map_;
	std::unordered_map<Partial, double> partial_to_first_derivative_sign_map_;
	std::unordered_map<Partial, double> partial_to_first_derivative_value_map_;
	std::unordered_map<std::pair<PartialCoordinateType, Partial>, double, Utils::PairHash, Utils::PairEquals> partial_to_second_derivative_value_map_;
};

#endif