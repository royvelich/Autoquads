#pragma once
#ifndef OPTIMIZATION_LIB_EDGE_PAIR_INTEGER_TRANSLATION_OBJECTIVE_H
#define OPTIMIZATION_LIB_EDGE_PAIR_INTEGER_TRANSLATION_OBJECTIVE_H

// Eigen includes
#include <Eigen/Core>

// Optimization lib includes
#include "../summation_objective.h"
#include "../coordinate_diff_objective.h"
#include "../cross_coordinate_diff_objective.h"
#include "../periodic_objective.h"
#include "../../data_providers/edge_pair_data_provider.h"
#include "../../data_providers/empty_data_provider.h"

template <Eigen::StorageOptions StorageOrder_>
class EdgePairIntegerTranslationObjective : public SummationObjective<PeriodicObjective<StorageOrder_>, Eigen::SparseVector<double>>
{
public:
	/**
	 * Constructors and destructor
	 */
	EdgePairIntegerTranslationObjective(const std::shared_ptr<MeshDataProvider>& mesh_data_provider, const std::shared_ptr<EdgePairDataProvider>& edge_pair_data_provider, const bool enforce_children_psd = true) :
		SummationObjective(mesh_data_provider, edge_pair_data_provider, "Edge Pair Translation Objective", enforce_children_psd)
	{
		this->Initialize();
	}

	virtual ~EdgePairIntegerTranslationObjective()
	{

	}

	/**
	 * Setters
	 */
	void SetInterval(const double interval)
	{
		for (auto& periodic_objective : periodic_objectives)
		{
			periodic_objective->SetPeriod(interval);
		}
	}

	/**
	 * Getters
	 */
	const EdgePairDataProvider& GetEdgePairDataProvider() const
	{
		return *std::dynamic_pointer_cast<EdgePairDataProvider>(this->data_provider_);
	}

	/**
	 * Public overrides
	 */
	void PreInitialize() override
	{
		auto edge_pair_data_provider = std::dynamic_pointer_cast<EdgePairDataProvider>(this->data_provider_);
		auto empty_data_provider = std::make_shared<EmptyDataProvider>(this->GetMeshDataProvider());

		auto v1_x_coordinate_diff_data_provider = std::make_shared<CoordinateDiffDataProvider>(this->mesh_data_provider_, edge_pair_data_provider->GetEdge1Vertex1Index(), edge_pair_data_provider->GetEdge2Vertex1Index(), RDS::CoordinateType::X);
		auto v1_y_coordinate_diff_data_provider = std::make_shared<CoordinateDiffDataProvider>(this->mesh_data_provider_, edge_pair_data_provider->GetEdge1Vertex1Index(), edge_pair_data_provider->GetEdge2Vertex1Index(), RDS::CoordinateType::Y);
		//auto v2_x_coordinate_diff_data_provider = std::make_shared<CoordinateDiffDataProvider>(this->mesh_data_provider_, edge_pair_data_provider->GetEdge1Vertex2Index(), edge_pair_data_provider->GetEdge2Vertex2Index(), RDS::CoordinateType::X);
		//auto v2_y_coordinate_diff_data_provider = std::make_shared<CoordinateDiffDataProvider>(this->mesh_data_provider_, edge_pair_data_provider->GetEdge1Vertex2Index(), edge_pair_data_provider->GetEdge2Vertex2Index(), RDS::CoordinateType::Y);

		auto v1_x_coordinate_diff_objective = std::make_shared<CoordinateDiffObjective<StorageOrder_>>(this->GetMeshDataProvider(), v1_x_coordinate_diff_data_provider);
		auto v1_y_coordinate_diff_objective = std::make_shared<CoordinateDiffObjective<StorageOrder_>>(this->GetMeshDataProvider(), v1_y_coordinate_diff_data_provider);
		//auto v2_x_coordinate_diff_objective = std::make_shared<CoordinateDiffObjective<StorageOrder_>>(this->GetMeshDataProvider(), v2_x_coordinate_diff_data_provider);
		//auto v2_y_coordinate_diff_objective = std::make_shared<CoordinateDiffObjective<StorageOrder_>>(this->GetMeshDataProvider(), v2_y_coordinate_diff_data_provider);

		auto periodic_v1_x_coordinate_diff_objective = std::make_shared<PeriodicObjective<StorageOrder_>>(this->GetMeshDataProvider(), empty_data_provider, v1_x_coordinate_diff_objective, 1.0f, this->GetEnforceChildrenPsd());
		auto periodic_v1_y_coordinate_diff_objective = std::make_shared<PeriodicObjective<StorageOrder_>>(this->GetMeshDataProvider(), empty_data_provider, v1_y_coordinate_diff_objective, 1.0f, this->GetEnforceChildrenPsd());
		//auto periodic_v2_x_coordinate_diff_objective = std::make_shared<PeriodicObjective<StorageOrder_>>(this->GetMeshDataProvider(), empty_data_provider, v2_x_coordinate_diff_objective, 1.0f, this->GetEnforceChildrenPsd());
		//auto periodic_v2_y_coordinate_diff_objective = std::make_shared<PeriodicObjective<StorageOrder_>>(this->GetMeshDataProvider(), empty_data_provider, v2_y_coordinate_diff_objective, 1.0f, this->GetEnforceChildrenPsd());

		this->AddObjectiveFunction(periodic_v1_x_coordinate_diff_objective);
		this->AddObjectiveFunction(periodic_v1_y_coordinate_diff_objective);
		//this->AddObjectiveFunction(periodic_v2_x_coordinate_diff_objective);
		//this->AddObjectiveFunction(periodic_v2_y_coordinate_diff_objective);

		periodic_objectives.push_back(periodic_v1_x_coordinate_diff_objective);
		periodic_objectives.push_back(periodic_v1_y_coordinate_diff_objective);
		//periodic_objectives.push_back(periodic_v2_x_coordinate_diff_objective);
		//periodic_objectives.push_back(periodic_v2_y_coordinate_diff_objective);
	}

private:
	/**
	 * Private fields
	 */
	std::vector<std::shared_ptr<PeriodicObjective<StorageOrder_>>> periodic_objectives;
};

#endif