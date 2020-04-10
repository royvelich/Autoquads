#pragma once
#ifndef OPTIMIZATION_LIB_SEAMLESS_OBJECTIVE_H
#define OPTIMIZATION_LIB_SEAMLESS_OBJECTIVE_H

// C includes
#define _USE_MATH_DEFINES
#include <math.h>

// STL includes
#include <vector>

// TBB includes
#include <tbb/concurrent_vector.h>

// Eigen includes
#include <Eigen/Core>

// Optimization lib includes
#include "../data_providers/empty_data_provider.h"
#include "../data_providers/edge_pair_data_provider.h"
#include "./summation_objective.h"
#include "./edge_pair/edge_pair_angle_objective.h"
#include "./edge_pair/edge_pair_length_objective.h"
#include "./edge_pair/edge_pair_integer_translation_objective.h"
#include "./edge_pair/edge_pair_translation_objective.h"
#include "../objective_functions/periodic_objective.h"

template <Eigen::StorageOptions StorageOrder_>
class SeamlessObjective : public SummationObjective<ObjectiveFunction<StorageOrder_, Eigen::SparseVector<double>>, Eigen::VectorXd>
{
public:
	/**
	 * Public type definitions
	 */
	enum class Properties : int32_t
	{
		Zeta = SummationObjective<SparseObjectiveFunction<StorageOrder_>, Eigen::VectorXd>::Properties::Count_,
		AngleValuePerEdge,
		LengthValuePerEdge,
		EdgeAngleWeight,
		EdgeLengthWeight,

		AngleWeight,
		LengthWeight,
		TranslationWeight,
		
		Interval
	};
	
	/**
	 * Constructors and destructor
	 */
	SeamlessObjective(const std::shared_ptr<MeshDataProvider>& mesh_data_provider, const std::shared_ptr<EmptyDataProvider>& empty_data_provider, const std::string& name, const bool enforce_children_psd = true) :
		SummationObjective(mesh_data_provider, empty_data_provider, name, enforce_children_psd),
		zeta_(1),
		angle_weight_(0),
		length_weight_(0),
		translation_weight_(0)
	{

	}

	SeamlessObjective(const std::shared_ptr<MeshDataProvider>& mesh_data_provider, const std::shared_ptr<EmptyDataProvider>& empty_data_provider, const bool enforce_children_psd = true) :
		SeamlessObjective(mesh_data_provider, empty_data_provider, "Seamless", enforce_children_psd)
	{

	}

	virtual ~SeamlessObjective()
	{

	}

	/**
	 * Setters
	 */
	void SetInterval(const double interval)
	{
		for (auto& edge_pair_integer_translation_objective : edge_pair_integer_translation_objectives)
		{
			edge_pair_integer_translation_objective->SetInterval(interval);
		}

		interval_ = interval;
	}
	
	void SetZeta(const double zeta)
	{
		//for(auto& periodic_edge_pair_angle_objective : periodic_edge_pair_angle_objectives)
		//{
		//	periodic_edge_pair_angle_objective->SetWeight(1.0 - zeta);
		//}

		//for (auto& edge_pair_length_objective : edge_pair_length_objectives)
		//{
		//	edge_pair_length_objective->SetWeight(zeta);
		//}
		//
		//zeta_ = zeta;
	}

	void SetEdgeAngleWeight(const RDS::EdgeIndex edge_index, const double weight)
	{
		for (auto& periodic_edge_pair_angle_objective : periodic_edge_pair_angle_objectives)
		{
			auto edge_pair_angle_objective = std::dynamic_pointer_cast<EdgePairAngleObjective<StorageOrder_>>(periodic_edge_pair_angle_objective->GetInnerObjective());
			if (edge_pair_angle_objective->GetEdgePairDataProvider().GetDomainEdgeIndex() == edge_index)
			{
				edge_pair_angle_objective->SetWeight(weight);
			}
		}
	}

	void SetEdgeLengthWeight(const RDS::EdgeIndex edge_index, const double weight)
	{
		for (auto& edge_pair_length_objective : edge_pair_length_objectives)
		{
			if (edge_pair_length_objective->GetEdgePairDataProvider().GetDomainEdgeIndex() == edge_index)
			{
				edge_pair_length_objective->SetWeight(weight);
			}
		}
	}
	
	void SetAngleWeight(const double weight)
	{
		angle_weight_ = weight;
		for (auto& periodic_edge_pair_angle_objective : periodic_edge_pair_angle_objectives)
		{
			periodic_edge_pair_angle_objective->SetWeight(weight);
		}
	}

	void SetLengthWeight(const double weight)
	{
		length_weight_ = weight;
		for (auto& edge_pair_length_objective : edge_pair_length_objectives)
		{
			edge_pair_length_objective->SetWeight(weight);
		}
	}

	void SetTranslationWeight(const double weight)
	{
		translation_weight_ = weight;
		for (auto& edge_pair_integer_translation_objective : edge_pair_integer_translation_objectives)
		{
			edge_pair_integer_translation_objective->SetWeight(weight);
		}

		//for (auto& edge_pair_translation_objective : edge_pair_translation_objectives)
		//{
		//	edge_pair_translation_objective->SetWeight(weight);
		//}
	}

	bool SetProperty(const int32_t property_id, const std::any property_context, const std::any property_value) override
	{
		if (SummationObjective<ObjectiveFunction<StorageOrder_, Eigen::SparseVector<double>>, Eigen::VectorXd>::SetProperty(property_id, property_context, property_value))
		{
			return true;
		}

		const Properties properties = static_cast<Properties>(property_id);
		switch (properties)
		{
		case Properties::Zeta:
			SetZeta(std::any_cast<const double>(property_value));
			return true;
		case Properties::EdgeAngleWeight:
			SetEdgeAngleWeight(static_cast<RDS::EdgeIndex>(std::any_cast<double>(property_context)), std::any_cast<const double>(property_value));
			return true;
		case Properties::EdgeLengthWeight:
			SetEdgeLengthWeight(static_cast<RDS::EdgeIndex>(std::any_cast<double>(property_context)), std::any_cast<const double>(property_value));
			return true;
		case Properties::AngleWeight:
			SetAngleWeight(std::any_cast<const double>(property_value));
			return true;
		case Properties::LengthWeight:
			SetLengthWeight(std::any_cast<const double>(property_value));
			return true;
		case Properties::TranslationWeight:
			SetTranslationWeight(std::any_cast<const double>(property_value));
			return true;	
		case Properties::Interval:
			SetInterval(std::any_cast<const double>(property_value));
			return true;
		}

		return false;
	}

	/**
	 * Getters
	 */
	double GetZeta() const
	{
		return zeta_;
	}

	const Eigen::VectorXd& SeamlessObjective::GetAngleValuePerEdge(const ObjectiveFunctionBase::PropertyModifiers property_modifiers) const
	{
		switch (property_modifiers)
		{
		case ObjectiveFunctionBase::PropertyModifiers::Domain:
			return GetDomainAngleValuePerEdge();
		case ObjectiveFunctionBase::PropertyModifiers::Image:
			return GetImageAngleValuePerEdge();
		}
	}

	const Eigen::VectorXd& SeamlessObjective::GetLengthValuePerEdge(const ObjectiveFunctionBase::PropertyModifiers property_modifiers) const
	{
		switch (property_modifiers)
		{
		case ObjectiveFunctionBase::PropertyModifiers::Domain:
			return GetDomainLengthValuePerEdge();
		case ObjectiveFunctionBase::PropertyModifiers::Image:
			return GetImageLengthValuePerEdge();
		}
	}

	const Eigen::VectorXd& SeamlessObjective::GetImageAngleValuePerEdge() const
	{
		return image_angle_value_per_edge_;
	}

	const Eigen::VectorXd& SeamlessObjective::GetImageLengthValuePerEdge() const
	{
		return image_length_value_per_edge_;
	}

	const Eigen::VectorXd& SeamlessObjective::GetDomainAngleValuePerEdge() const
	{
		return domain_angle_value_per_edge_;
	}

	const Eigen::VectorXd& SeamlessObjective::GetDomainLengthValuePerEdge() const
	{
		return domain_length_value_per_edge_;
	}

	double GetEdgeAngleWeight(const RDS::EdgeIndex edge_index)
	{
		for (auto& periodic_edge_pair_angle_objective : periodic_edge_pair_angle_objectives)
		{
			auto edge_pair_angle_objective = std::dynamic_pointer_cast<EdgePairAngleObjective<StorageOrder_>>(periodic_edge_pair_angle_objective->GetInnerObjective());
			if (edge_pair_angle_objective->GetEdgePairDataProvider().GetDomainEdgeIndex() == edge_index)
			{
				return edge_pair_angle_objective->GetWeight();
			}
		}

		return 0;
	}

	double GetEdgeLengthWeight(const RDS::EdgeIndex edge_index)
	{
		for (auto& edge_pair_length_objective : edge_pair_length_objectives)
		{
			if (edge_pair_length_objective->GetEdgePairDataProvider().GetDomainEdgeIndex() == edge_index)
			{
				return edge_pair_length_objective->GetWeight();
			}
		}

		return 0;
	}

	double GetAngleWeight() const
	{
		return angle_weight_;
	}

	double GetLengthWeight() const
	{
		return length_weight_;
	}

	double GetTranslationWeight() const
	{
		return translation_weight_;
	}

	bool GetProperty(const int32_t property_id, const int32_t property_modifier_id, const std::any property_context, std::any& property_value) override
	{
		if (SummationObjective<ObjectiveFunction<StorageOrder_, Eigen::SparseVector<double>>, Eigen::VectorXd>::GetProperty(property_id, property_modifier_id, property_context, property_value))
		{
			return true;
		}

		const ObjectiveFunctionBase::PropertyModifiers property_modifiers = static_cast<ObjectiveFunctionBase::PropertyModifiers>(property_modifier_id);
		const Properties properties = static_cast<Properties>(property_id);
		switch (properties)
		{
		case Properties::Zeta:
			property_value = GetZeta();
			return true;
		case Properties::AngleValuePerEdge:
			property_value = GetAngleValuePerEdge(property_modifiers);
			return true;
		case Properties::LengthValuePerEdge:
			property_value = GetLengthValuePerEdge(property_modifiers);
			return true;
		case Properties::EdgeAngleWeight:
			property_value = GetEdgeAngleWeight(static_cast<RDS::EdgeIndex>(std::any_cast<double>(property_context)));
			return true;
		case Properties::EdgeLengthWeight:
			property_value = GetEdgeLengthWeight(static_cast<RDS::EdgeIndex>(std::any_cast<double>(property_context)));
			return true;
		case Properties::AngleWeight:
			property_value = GetAngleWeight();
			return true;
		case Properties::LengthWeight:
			property_value = GetLengthWeight();
			return true;
		case Properties::TranslationWeight:
			property_value = GetTranslationWeight();
			return true;		
		}

		return false;
	}

	/**
	 * Public methods
	 */	
	void AddEdgePairObjectives(const std::shared_ptr<EdgePairDataProvider>& edge_pair_data_provider)
	{	
		auto edge_pair_angle_objective = std::make_shared<EdgePairAngleObjective<StorageOrder_>>(this->GetMeshDataProvider(), edge_pair_data_provider, false);
		auto edge_pair_length_objective = std::make_shared<EdgePairLengthObjective<StorageOrder_>>(this->GetMeshDataProvider(), edge_pair_data_provider, false);
		auto edge_pair_integer_translation_objective = std::make_shared<EdgePairIntegerTranslationObjective<StorageOrder_>>(this->GetMeshDataProvider(), edge_pair_data_provider, this->GetEnforceChildrenPsd());
		//auto edge_pair_translation_objective = std::make_shared<EdgePairTranslationObjective<StorageOrder_>>(this->GetMeshDataProvider(), edge_pair_data_provider, this->GetEnforceChildrenPsd());
		
		double period = M_PI / 2;
		auto empty_data_provider = std::make_shared<EmptyDataProvider>(this->GetMeshDataProvider());
		std::shared_ptr<PeriodicObjective<StorageOrder_>> periodic_edge_pair_angle_objective = std::make_shared<PeriodicObjective<StorageOrder_>>(this->GetMeshDataProvider(), empty_data_provider, edge_pair_angle_objective, period, this->GetEnforceChildrenPsd());

		periodic_edge_pair_angle_objective->SetWeight(angle_weight_);
		edge_pair_length_objective->SetWeight(length_weight_);
		edge_pair_integer_translation_objective->SetWeight(translation_weight_);
		//edge_pair_translation_objective->SetWeight(translation_weight_);
		
		this->AddObjectiveFunction(periodic_edge_pair_angle_objective);
		this->AddObjectiveFunction(edge_pair_length_objective);
		this->AddObjectiveFunction(edge_pair_integer_translation_objective);
		//this->AddObjectiveFunction(edge_pair_translation_objective);
		
		periodic_edge_pair_angle_objectives.push_back(periodic_edge_pair_angle_objective);
		edge_pair_length_objectives.push_back(edge_pair_length_objective);
		edge_pair_angle_objectives.push_back(edge_pair_angle_objective);
		edge_pair_integer_translation_objectives.push_back(edge_pair_integer_translation_objective);
		//edge_pair_translation_objectives.push_back(edge_pair_translation_objective);
	}

protected:
	/**
	 * Protected overrides
	 */
	void PostInitialize() override
	{
		SummationObjective<ObjectiveFunction<StorageOrder_, Eigen::SparseVector<double>>, Eigen::VectorXd>::PostInitialize();
		image_angle_value_per_edge_.resize(this->mesh_data_provider_->GetImageEdgesCount());
		image_length_value_per_edge_.resize(this->mesh_data_provider_->GetImageEdgesCount());
		domain_angle_value_per_edge_.resize(this->mesh_data_provider_->GetDomainEdgesCount());
		domain_length_value_per_edge_.resize(this->mesh_data_provider_->GetDomainEdgesCount());
	}
	
private:
	
	/**
	 * Private overrides
	 */
	void CalculateValuePerEdge(Eigen::VectorXd& domain_value_per_edge, Eigen::VectorXd& image_value_per_edge) override
	{
		CalculateAngleValuePerEdge(domain_angle_value_per_edge_, image_angle_value_per_edge_);
		CalculateLengthValuePerEdge(domain_length_value_per_edge_, image_length_value_per_edge_);
 		domain_value_per_edge = domain_angle_value_per_edge_ + domain_length_value_per_edge_;
		image_value_per_edge = image_angle_value_per_edge_ + image_length_value_per_edge_;
	}
	
	/**
	 * Private methods
	 */
	void CalculateAngleValuePerEdge(Eigen::VectorXd& domain_angle_value_per_edge, Eigen::VectorXd& image_angle_value_per_edge)
	{
		domain_angle_value_per_edge.setZero();
		image_angle_value_per_edge.setZero();
		for(const auto& periodic_edge_pair_angle_objective : periodic_edge_pair_angle_objectives)
		{
			auto edge_pair_angle_objective = std::dynamic_pointer_cast<EdgePairAngleObjective<StorageOrder_>>(periodic_edge_pair_angle_objective->GetInnerObjective());
			const RDS::EdgeIndex domain_edge_index = edge_pair_angle_objective->GetEdgePairDataProvider().GetDomainEdgeIndex();
			const RDS::EdgeIndex image_edge_1_index = edge_pair_angle_objective->GetEdgePairDataProvider().GetImageEdge1Index();
			const RDS::EdgeIndex image_edge_2_index = edge_pair_angle_objective->GetEdgePairDataProvider().GetImageEdge2Index();

			const double value = periodic_edge_pair_angle_objective->GetValue();
			domain_angle_value_per_edge.coeffRef(domain_edge_index) += value;
			image_angle_value_per_edge.coeffRef(image_edge_1_index) += value;
			image_angle_value_per_edge.coeffRef(image_edge_2_index) += value;
		}
	}

	void CalculateLengthValuePerEdge(Eigen::VectorXd& domain_length_value_per_edge, Eigen::VectorXd& image_length_value_per_edge)
	{
		domain_length_value_per_edge.setZero();
		image_length_value_per_edge.setZero();
		for (const auto& edge_pair_length_objective : edge_pair_length_objectives)
		{
			const RDS::EdgeIndex domain_edge_index = edge_pair_length_objective->GetEdgePairDataProvider().GetDomainEdgeIndex();
			const RDS::EdgeIndex image_edge_1_index = edge_pair_length_objective->GetEdgePairDataProvider().GetImageEdge1Index();
			const RDS::EdgeIndex image_edge_2_index = edge_pair_length_objective->GetEdgePairDataProvider().GetImageEdge2Index();

			const double value = edge_pair_length_objective->GetValue();
			domain_length_value_per_edge.coeffRef(domain_edge_index) += value;
			image_length_value_per_edge.coeffRef(image_edge_1_index) += value;
			image_length_value_per_edge.coeffRef(image_edge_2_index) += value;
		}
	}
	
	/**
	 * Private fields
	 */
	double zeta_;
	double interval_;

	double angle_weight_;
	double length_weight_;
	double translation_weight_;
	
	tbb::concurrent_vector<std::shared_ptr<EdgePairLengthObjective<StorageOrder_>>> edge_pair_length_objectives;
	tbb::concurrent_vector<std::shared_ptr<EdgePairAngleObjective<StorageOrder_>>> edge_pair_angle_objectives;
	tbb::concurrent_vector<std::shared_ptr<PeriodicObjective<StorageOrder_>>> periodic_edge_pair_angle_objectives;
	tbb::concurrent_vector<std::shared_ptr<EdgePairIntegerTranslationObjective<StorageOrder_>>> edge_pair_integer_translation_objectives;
	tbb::concurrent_vector<std::shared_ptr<EdgePairTranslationObjective<StorageOrder_>>> edge_pair_translation_objectives;
	
	Eigen::VectorXd image_angle_value_per_edge_;
	Eigen::VectorXd image_length_value_per_edge_;
	Eigen::VectorXd domain_angle_value_per_edge_;
	Eigen::VectorXd domain_length_value_per_edge_;
};

#endif