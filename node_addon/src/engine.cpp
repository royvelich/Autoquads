// STL includes
#include <sstream>
#include <any>

// LIBIGL includes
#include <igl/readOFF.h>
#include <igl/readOBJ.h>

// Optimization lib includes
#include "../include/engine.h"
#include "libs/optimization_lib/include/objective_functions/objective_function.h"
#include "libs/optimization_lib/include/objective_functions/dense_objective_function.h"
#include "libs/optimization_lib/include/objective_functions/position/face_barycenter_position_objective.h"
#include "libs/optimization_lib/include/objective_functions/position/face_vertices_position_objective.h"
#include "libs/optimization_lib/include/objective_functions/edge_pair/edge_pair_angle_objective.h"
#include "libs/optimization_lib/include/objective_functions/edge_pair/edge_pair_length_objective.h"
#include "libs/optimization_lib/include/objective_functions/singularity/singular_points_position_objective.h"

Napi::FunctionReference Engine::constructor;

Napi::Object Engine::Init(Napi::Env env, Napi::Object exports)
{
	Napi::HandleScope scope(env);

	Napi::Function func = DefineClass(env, "Engine", {
		InstanceMethod("loadModel", &Engine::LoadModel),
		InstanceMethod("resumeSolver", &Engine::ResumeSolver),
		InstanceMethod("pauseSolver", &Engine::PauseSolver),
		InstanceMethod("constrainFacePosition", &Engine::ConstrainFacePosition),
		InstanceMethod("updateConstrainedFacePosition", &Engine::UpdateConstrainedFacePosition),
		InstanceMethod("unconstrainFacePosition", &Engine::UnconstrainFacePosition),
		InstanceMethod("getDomainFacesCount", &Engine::GetDomainFacesCount),
		InstanceMethod("getImageFacesCount", &Engine::GetImageFacesCount),
		InstanceMethod("getDomainEdgesCount", &Engine::GetDomainEdgesCount),
		InstanceMethod("getImageEdgesCount", &Engine::GetImageEdgesCount),
		InstanceMethod("getDomainVerticesCount", &Engine::GetDomainVerticesCount),
		InstanceMethod("getImageVerticesCount", &Engine::GetImageVerticesCount),
		InstanceMethod("getDomainFaces", &Engine::GetDomainFaces),
		InstanceMethod("getImageFaces", &Engine::GetImageFaces),
		InstanceMethod("getDomainEdges", &Engine::GetDomainEdges),
		InstanceMethod("getImageEdges", &Engine::GetImageEdges),
		InstanceMethod("getDomainVertices", &Engine::GetDomainVertices),
		InstanceMethod("getImageVertices", &Engine::GetImageVertices),
		InstanceMethod("getDomainBufferedFaces", &Engine::GetDomainBufferedFaces),
		InstanceMethod("getImageBufferedFaces", &Engine::GetImageBufferedFaces),
		InstanceMethod("getDomainBufferedEdges", &Engine::GetDomainBufferedEdges),
		InstanceMethod("getImageBufferedEdges", &Engine::GetImageBufferedEdges),
		InstanceMethod("getDomainBufferedVertices", &Engine::GetDomainBufferedVertices),
		InstanceMethod("getImageBufferedVertices", &Engine::GetImageBufferedVertices),
		InstanceMethod("getDomainBufferedUvs", &Engine::GetDomainBufferedUvs),
		InstanceMethod("getImageBufferedUvs", &Engine::GetImageBufferedUvs),
		InstanceMethod("getDomainFaceEdgeAdjacency", &Engine::GetDomainFaceEdgeAdjacency),
		InstanceMethod("getDomainEdgeFaceAdjacency", &Engine::GetDomainEdgeFaceAdjacency),
		InstanceMethod("getImageFaceEdgeAdjacency", &Engine::GetImageFaceEdgeAdjacency),
		InstanceMethod("getImageEdgeFaceAdjacency", &Engine::GetImageEdgeFaceAdjacency),	
		InstanceMethod("getObjectiveFunctionProperty", &Engine::GetObjectiveFunctionProperty),
		InstanceMethod("setObjectiveFunctionProperty", &Engine::SetObjectiveFunctionProperty),
		InstanceMethod("setAlgorithmType", &Engine::SetAlgorithmType),
		InstanceAccessor("positionWeight", &Engine::GetPositionWeight, &Engine::SetPositionWeight),
		InstanceAccessor("seamlessWeight", &Engine::GetSeamlessWeight, &Engine::SetSeamlessWeight),
		InstanceAccessor("lambda", &Engine::GetLambda, &Engine::SetLambda),
		InstanceAccessor("delta", &Engine::GetDelta, &Engine::SetDelta),
		InstanceAccessor("objectiveFunctionsData", &Engine::GetObjectiveFunctionsData, nullptr)
	});

	constructor = Napi::Persistent(func);
	constructor.SuppressDestruct();

	exports.Set("Engine", func);
	return exports;
}

Engine::Engine(const Napi::CallbackInfo& info) : 
	Napi::ObjectWrap<Engine>(info),
	mesh_wrapper_(std::make_shared<MeshWrapper>())
{
	properties_map_.insert({ "value", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::Value) });
	properties_map_.insert({ "value_per_vertex", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::ValuePerVertex) });
	properties_map_.insert({ "value_per_edge", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::ValuePerEdge) });
	properties_map_.insert({ "gradient", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::Gradient) });
	properties_map_.insert({ "gradient_norm", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::GradientNorm) });
	properties_map_.insert({ "hessian", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::Hessian) });
	properties_map_.insert({ "weight", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::Weight) });
	properties_map_.insert({ "name", static_cast<uint32_t>(ObjectiveFunctionBase::Properties::Name) });
	properties_map_.insert({ "delta", static_cast<uint32_t>(Separation<Eigen::StorageOptions::RowMajor>::Properties::Delta) });
	properties_map_.insert({ "zeta", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::Zeta) });
	properties_map_.insert({ "angle_value_per_edge", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::AngleValuePerEdge) });
	properties_map_.insert({ "length_value_per_edge", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::LengthValuePerEdge) });
	properties_map_.insert({ "edge_angle_weight", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::EdgeAngleWeight) });
	properties_map_.insert({ "edge_length_weight", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::EdgeLengthWeight) });
	properties_map_.insert({ "angle_weight", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::AngleWeight) });
	properties_map_.insert({ "length_weight", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::LengthWeight) });
	properties_map_.insert({ "translation_weight", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::TranslationWeight) });
	properties_map_.insert({ "translation_interval", static_cast<uint32_t>(SeamlessObjective<Eigen::StorageOptions::RowMajor>::Properties::Interval) });
	properties_map_.insert({ "interval", static_cast<uint32_t>(SingularPointsPositionObjective<Eigen::StorageOptions::RowMajor>::Properties::Interval) });
	properties_map_.insert({ "singularity_weight_per_vertex", static_cast<uint32_t>(SingularPointsPositionObjective<Eigen::StorageOptions::RowMajor>::Properties::SingularityWeightPerVertex) });
	properties_map_.insert({ "negative_angular_defect_singularities_indices", static_cast<uint32_t>(SingularPointsPositionObjective<Eigen::StorageOptions::RowMajor>::Properties::NegativeAngularDefectSingularitiesIndices) });
	properties_map_.insert({ "positive_angular_defect_singularities_indices", static_cast<uint32_t>(SingularPointsPositionObjective<Eigen::StorageOptions::RowMajor>::Properties::PositiveAngularDefectSingularitiesIndices) });

	property_modifiers_map_.insert({ "none", static_cast<uint32_t>(ObjectiveFunctionBase::PropertyModifiers::None) });
	property_modifiers_map_.insert({ "domain", static_cast<uint32_t>(ObjectiveFunctionBase::PropertyModifiers::Domain) });
	property_modifiers_map_.insert({ "image", static_cast<uint32_t>(ObjectiveFunctionBase::PropertyModifiers::Image) });
	
	empty_data_provider_ = std::make_shared<EmptyDataProvider>(mesh_wrapper_);
	plain_data_provider_ = std::make_shared<PlainDataProvider>(mesh_wrapper_);
	
	// TODO: Expose interface for addition and removal of objective function
	separation_ = std::make_shared<Separation<Eigen::StorageOptions::RowMajor>>(mesh_wrapper_, empty_data_provider_);
	symmetric_dirichlet_ = std::make_shared<SymmetricDirichlet<Eigen::StorageOptions::RowMajor>>(mesh_wrapper_, empty_data_provider_);
	seamless_ = std::make_shared<SeamlessObjective<Eigen::StorageOptions::RowMajor>>(mesh_wrapper_, empty_data_provider_);
	singular_points_ = std::make_shared<SingularPointsPositionObjective<Eigen::StorageOptions::RowMajor>>(mesh_wrapper_, empty_data_provider_, 1);
  	position_ = std::make_shared<SummationObjective<ObjectiveFunction<Eigen::StorageOptions::RowMajor, Eigen::VectorXd>, Eigen::VectorXd>>(mesh_wrapper_, empty_data_provider_, std::string("Position"));

	objective_functions_.push_back(separation_);
	objective_functions_.push_back(symmetric_dirichlet_);
	objective_functions_.push_back(seamless_);
	objective_functions_.push_back(singular_points_);
	objective_functions_.push_back(position_);
	
	autoquads_objective_functions_.push_back(separation_);
	autoquads_objective_functions_.push_back(symmetric_dirichlet_);
	autoquads_objective_functions_.push_back(seamless_);
	autoquads_objective_functions_.push_back(singular_points_);

	//autocuts_objective_functions_.push_back(separation_);
	autocuts_objective_functions_.push_back(symmetric_dirichlet_);
	autocuts_objective_functions_.push_back(seamless_);
	autocuts_objective_functions_.push_back(singular_points_);
	autocuts_objective_functions_.push_back(position_);
	
	autocuts_summation_objective_ = std::make_shared<SummationObjective<ObjectiveFunction<Eigen::StorageOptions::RowMajor, Eigen::VectorXd>, Eigen::VectorXd>>(mesh_wrapper_, empty_data_provider_, autocuts_objective_functions_, false);
	autoquads_summation_objective_ = std::make_shared<SummationObjective<ObjectiveFunction<Eigen::StorageOptions::RowMajor, Eigen::VectorXd>, Eigen::VectorXd>>(mesh_wrapper_, empty_data_provider_, autoquads_objective_functions_, false);
	summation_objectives_.push_back(autocuts_summation_objective_);
	summation_objectives_.push_back(autoquads_summation_objective_);
	
	summation_objective_ = autocuts_summation_objective_;
	
	mesh_wrapper_->RegisterModelLoadedCallback([this]() {
		/**
		 * Initialize objective functions
		 */
		face_data_providers_.clear();
		face_to_face_data_provider_map_.clear();
		for (auto& face : mesh_wrapper_->GetImageFacesSTL())
		{
			auto face_data_provider = std::make_shared<FaceDataProvider>(mesh_wrapper_, face);
			face_data_providers_.push_back(face_data_provider);
			face_to_face_data_provider_map_.insert(std::make_pair(face, face_data_provider));
		}

		const auto edge_pair_descriptors = mesh_wrapper_->GetEdgePairDescriptors();
		edge_pair_data_providers_.resize(edge_pair_descriptors.size());
		#pragma omp parallel for
		for (int64_t i = 0; i < edge_pair_descriptors.size(); i++)
		{
			auto edge_pair_descriptor = edge_pair_descriptors[i];
			const auto edge_pair_data_provider = std::make_shared<EdgePairDataProvider>(mesh_wrapper_, edge_pair_descriptor);
			edge_pair_data_providers_[i] = edge_pair_data_provider;
		}

		const auto face_fans = mesh_wrapper_->GetFaceFans();
		face_fan_data_providers_.resize(face_fans.size());
		#pragma omp parallel for
		for (int64_t i = 0; i < face_fans.size(); i++)
		{
			auto face_fan = face_fans[i];
			const auto face_fan_data_provider = std::make_shared<FaceFanDataProvider>(mesh_wrapper_, face_fan);
			face_fan_data_providers_[i] = face_fan_data_provider;
		}

		#pragma omp parallel for
		for (int64_t i = 0; i < edge_pair_data_providers_.size(); i++)
		{
			seamless_->AddEdgePairObjectives(edge_pair_data_providers_[i]);
		}

		#pragma omp parallel for
		for (int64_t i = 0; i < face_fan_data_providers_.size(); i++)
		{
			singular_points_->AddSingularPointObjective(face_fan_data_providers_[i]);
		}

		autocuts_summation_objective_->Initialize();
		//autoquads_summation_objective_->Initialize();
		
		/**
		 * Create newton method iterator
		 */
		auto image_vertices = mesh_wrapper_->GetImageVertices();
		auto x0 = Eigen::Map<const Eigen::VectorXd>(image_vertices.data(), image_vertices.cols() * image_vertices.rows());
		newton_method_ = std::make_unique<NewtonMethod<PardisoSolver, Eigen::StorageOptions::RowMajor>>(summation_objective_, x0);
		newton_method_->EnableFlipAvoidingLineSearch(mesh_wrapper_->GetImageFaces());
	});
}


Napi::Value Engine::GetDomainVerticesCount(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number count = Napi::Number::New(env, mesh_wrapper_->GetDomainVertices().rows());

	return count;
}

Napi::Value Engine::GetImageVerticesCount(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number count = Napi::Number::New(env, mesh_wrapper_->GetImageVertices().rows());

	return count;
}

Napi::Value Engine::GetDomainFacesCount(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number count = Napi::Number::New(env, mesh_wrapper_->GetDomainFaces().rows());

	return count;
}

Napi::Value Engine::GetImageFacesCount(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number count = Napi::Number::New(env, mesh_wrapper_->GetImageFaces().rows());

	return count;
}

Napi::Value Engine::GetDomainEdgesCount(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number count = Napi::Number::New(env, mesh_wrapper_->GetDomainEdges().rows());

	return count;
}

Napi::Value Engine::GetImageEdgesCount(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number count = Napi::Number::New(env, mesh_wrapper_->GetImageEdges().rows());

	return count;
}

Napi::Value Engine::GetDomainFaces(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	return CreateFaces(env, mesh_wrapper_->GetDomainFaces());
}

Napi::Value Engine::GetImageFaces(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	return CreateFaces(env, mesh_wrapper_->GetImageFaces());
}

Napi::Value Engine::GetDomainEdges(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	return CreateEdges(env, mesh_wrapper_->GetDomainEdges());
}

Napi::Value Engine::GetImageEdges(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	return CreateEdges(env, mesh_wrapper_->GetImageEdges());
}

Napi::Value Engine::GetDomainVertices(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	return CreateVerticesArray(env, mesh_wrapper_->GetDomainVertices());
}

Napi::Value Engine::GetImageVertices(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	TryUpdateImageVertices();

	return CreateVerticesArray(env, mesh_wrapper_->GetImageVertices());
}

Napi::Value Engine::GetDomainBufferedFaces(const Napi::CallbackInfo& info)
{
	return GetBufferedFaces(info, FacesSource::DOMAIN_FACES);
}

Napi::Value Engine::GetImageBufferedFaces(const Napi::CallbackInfo& info)
{
	return GetBufferedFaces(info, FacesSource::IMAGE_FACES);
}

Napi::Value Engine::GetDomainBufferedEdges(const Napi::CallbackInfo& info)
{
	return GetBufferedEdges(info, EdgesSource::DOMAIN_EDGES);
}

Napi::Value Engine::GetImageBufferedEdges(const Napi::CallbackInfo& info)
{
	return GetBufferedEdges(info, EdgesSource::IMAGE_EDGES);
}

Napi::Value Engine::GetDomainBufferedVertices(const Napi::CallbackInfo& info)
{
	TryUpdateImageVertices();
	return GetBufferedVertices(info, VerticesSource::DOMAIN_VERTICES);
}

Napi::Value Engine::GetImageBufferedVertices(const Napi::CallbackInfo& info)
{
	TryUpdateImageVertices();
	return GetBufferedVertices(info, VerticesSource::IMAGE_VERTICES);
}

Napi::Value Engine::GetDomainBufferedUvs(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	TryUpdateImageVertices();

	return CreateBufferedUvsArray(env, mesh_wrapper_->GetImageVertices(), mesh_wrapper_->GetImageFaces());
}

Napi::Value Engine::GetImageBufferedUvs(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	TryUpdateImageVertices();

	return CreateBufferedUvsArray(env, mesh_wrapper_->GetImageVertices(), mesh_wrapper_->GetImageFaces());
}

Napi::Int32Array Engine::GetBufferedFaces(const Napi::CallbackInfo& info, const FacesSource faces_source) const
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Create buffered faces array
	 */
	switch (faces_source)
	{
	case FacesSource::DOMAIN_FACES:
		return CreateBufferedFacesArray(env, mesh_wrapper_->GetDomainFaces());
	case FacesSource::IMAGE_FACES:
		return CreateBufferedFacesArray(env, mesh_wrapper_->GetImageFaces());
	}

	Napi::TypeError::New(env, "Unknown buffered primitive type").ThrowAsJavaScriptException();
	return Napi::Int32Array::New(env, 0);
}

Napi::Int32Array Engine::GetBufferedEdges(const Napi::CallbackInfo& info, const EdgesSource edges_source) const
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Create buffered edges array
	 */
	switch (edges_source)
	{
	case EdgesSource::DOMAIN_EDGES:
		return CreateBufferedEdgesArray(env, mesh_wrapper_->GetDomainEdges());
	case EdgesSource::IMAGE_EDGES:
		return CreateBufferedEdgesArray(env, mesh_wrapper_->GetImageEdges());
	}

	Napi::TypeError::New(env, "Unknown buffered primitive type").ThrowAsJavaScriptException();
	return Napi::Int32Array::New(env, 0);
}

Napi::Float32Array Engine::GetBufferedVertices(const Napi::CallbackInfo& info, const VerticesSource vertices_source)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (info.Length() >= 1)
	{
		if (!info[0].IsNumber())
		{
			Napi::TypeError::New(env, "First argument is expected to be a Number").ThrowAsJavaScriptException();
			return Napi::Float32Array::New(env, 0);
		}
	}
	else
	{
		Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
		return Napi::Float32Array::New(env, 0);
	}

	/**
	 * Create buffered vertices array according to the provided primitive type
	 */
	Napi::Number argument1 = info[0].As<Napi::Number>();
	BufferedPrimitiveType buffered_primitive_type = static_cast<BufferedPrimitiveType>(argument1.Uint32Value());

	switch (buffered_primitive_type)
	{
	case BufferedPrimitiveType::TRIANGLE:
		switch (vertices_source)
		{
		case VerticesSource::DOMAIN_VERTICES:
			return CreateBufferedVerticesArray(env, mesh_wrapper_->GetDomainVertices(), mesh_wrapper_->GetDomainFaces());
		case VerticesSource::IMAGE_VERTICES:
			return CreateBufferedVerticesArray(env, mesh_wrapper_->GetImageVertices(), mesh_wrapper_->GetImageFaces());
		}
		break;

	case BufferedPrimitiveType::EDGE:
		switch (vertices_source)
		{
		case VerticesSource::DOMAIN_VERTICES:
			return CreateBufferedVerticesArray(env, mesh_wrapper_->GetDomainVertices(), mesh_wrapper_->GetDomainEdges());
		case VerticesSource::IMAGE_VERTICES:
			return CreateBufferedVerticesArray(env, mesh_wrapper_->GetImageVertices(), mesh_wrapper_->GetImageEdges());
		}
		break;
		
	case BufferedPrimitiveType::VERTEX:
		switch (vertices_source)
		{
		case VerticesSource::DOMAIN_VERTICES:
			return CreateBufferedVerticesArray(env, mesh_wrapper_->GetDomainVertices());
		case VerticesSource::IMAGE_VERTICES:
			return CreateBufferedVerticesArray(env, mesh_wrapper_->GetImageVertices());
		}
		break;
	}

	Napi::TypeError::New(env, "Unknown buffered primitive type").ThrowAsJavaScriptException();
	return Napi::Float32Array::New(env, 0);
}

Napi::Value Engine::GetFaceEdgeAdjacency(const Napi::CallbackInfo& info, const DataSource data_source)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);
	
	MeshWrapper::FI2EIsMap fi_2_ei;
	switch(data_source)
	{
	case DataSource::DOMAIN_DATA:
		fi_2_ei = mesh_wrapper_->GetDomainFaceEdgeAdjacency();
		break;
	case DataSource::IMAGE_DATA:
		fi_2_ei = mesh_wrapper_->GetImageFaceEdgeAdjacency();
		break;
	}

	Napi::Object face_edge_adjacency_object = Napi::Object::New(env);
	for (const auto& adjacency_entry : fi_2_ei)
	{
		Napi::Array adjacency_list_array = Napi::Array::New(env, adjacency_entry.second.size());
		face_edge_adjacency_object.Set(adjacency_entry.first, adjacency_list_array);

		for(std::size_t i = 0; i < adjacency_entry.second.size(); i++)
		{
			adjacency_list_array[i] = Napi::Number::New(env, adjacency_entry.second[i]);
		}
	}

	return face_edge_adjacency_object;
}

Napi::Value Engine::GetEdgeFaceAdjacency(const Napi::CallbackInfo& info, const DataSource data_source)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	MeshWrapper::EI2FIsMap ei_2_fi;
	switch (data_source)
	{
	case DataSource::DOMAIN_DATA:
		ei_2_fi = mesh_wrapper_->GetDomainEdgeFaceAdjacency();
		break;
	case DataSource::IMAGE_DATA:
		ei_2_fi = mesh_wrapper_->GetImageEdgeFaceAdjacency();
		break;
	}

	Napi::Object edge_face_adjacency_object = Napi::Object::New(env);
	for (const auto& adjacency_entry : ei_2_fi)
	{
		Napi::Array adjacency_list_array = Napi::Array::New(env, adjacency_entry.second.size());
		edge_face_adjacency_object.Set(adjacency_entry.first, adjacency_list_array);

		for (std::size_t i = 0; i < adjacency_entry.second.size(); i++)
		{
			adjacency_list_array[i] = Napi::Number::New(env, adjacency_entry.second[i]);
		}
	}

	return edge_face_adjacency_object;
}

Napi::Value Engine::GetDomainFaceEdgeAdjacency(const Napi::CallbackInfo& info)
{
	return GetFaceEdgeAdjacency(info, DataSource::DOMAIN_DATA);
}

Napi::Value Engine::GetDomainEdgeFaceAdjacency(const Napi::CallbackInfo& info)
{
	return GetEdgeFaceAdjacency(info, DataSource::DOMAIN_DATA);
}

Napi::Value Engine::GetImageFaceEdgeAdjacency(const Napi::CallbackInfo& info)
{
	return GetFaceEdgeAdjacency(info, DataSource::IMAGE_DATA);
}

Napi::Value Engine::GetImageEdgeFaceAdjacency(const Napi::CallbackInfo& info)
{
	return GetEdgeFaceAdjacency(info, DataSource::IMAGE_DATA);
}

std::shared_ptr<ObjectiveFunction<Eigen::StorageOptions::RowMajor, Eigen::VectorXd>> Engine::GetObjectiveFunctionByName(const std::string& name)
{
	for(const auto& summation_objective : summation_objectives_)
	{
		auto objective_function = summation_objective->GetObjectiveFunction(name);
		if(objective_function != nullptr)
		{
			return objective_function;
		}
	}

	return nullptr;
}

Napi::Value Engine::GetObjectiveFunctionProperty(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (info.Length() >= 4)
	{
		if (!info[0].IsString())
		{
			Napi::TypeError::New(env, "First argument is expected to be a String").ThrowAsJavaScriptException();
			return Napi::Value();
		}

		if (!info[1].IsString())
		{
			Napi::TypeError::New(env, "Second argument is expected to be a String").ThrowAsJavaScriptException();
			return Napi::Value();
		}

		if (!info[2].IsString())
		{
			Napi::TypeError::New(env, "Third argument is expected to be a String").ThrowAsJavaScriptException();
			return Napi::Value();
		}
	}
	else
	{
		Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Get objective function by name
	 */
	const std::string objective_function_name = info[0].ToString();
	const auto objective_function = GetObjectiveFunctionByName(objective_function_name);
	if(objective_function == nullptr)
	{
		Napi::TypeError::New(env, "Objective function could not be found").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Get property
	 */
	const std::string property_name = info[1].ToString();
	const std::string property_modifier_name = info[2].ToString();
	if(properties_map_.find(property_name) != properties_map_.end())
	{
		const uint32_t property_id = properties_map_.at(property_name);
		const uint32_t property_modifier_id = property_modifiers_map_.at(property_modifier_name);
		std::any any_value;
		if (objective_function->GetProperty(property_id, property_modifier_id, JSToNative(env, info[3]), any_value))
		{
			return NativeToJS(env, any_value);
		}
		
		Napi::TypeError::New(env, "Couldn't get property").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Property not found
	 */
	Napi::TypeError::New(env, "Property name could not be found").ThrowAsJavaScriptException();
	return Napi::Value();
}

Napi::Value Engine::SetObjectiveFunctionProperty(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (info.Length() >= 3)
	{
		if (!info[0].IsString())
		{
			Napi::TypeError::New(env, "First argument is expected to be a String").ThrowAsJavaScriptException();
			return Napi::Value();
		}

		if (!info[1].IsString())
		{
			Napi::TypeError::New(env, "Second argument is expected to be a String").ThrowAsJavaScriptException();
			return Napi::Value();
		}
	}
	else
	{
		Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Get objective function by name
	 */
	const std::string objective_function_name = info[0].ToString();
	const auto objective_function = GetObjectiveFunctionByName(objective_function_name);
	if (objective_function == nullptr)
	{
		Napi::TypeError::New(env, "Objective function could not be found").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Set property
	 */
	const std::string property_name = info[1].ToString();
	if (properties_map_.find(property_name) != properties_map_.end())
	{
		const uint32_t property_id = properties_map_.at(property_name);
		if(!objective_function->SetProperty(property_id, JSToNative(env, info[2]), JSToNative(env, info[3])))
		{
			Napi::TypeError::New(env, "Couldn't set property").ThrowAsJavaScriptException();
		}

		return Napi::Value();
	}

	/**
	 * Property not found
	 */
	Napi::TypeError::New(env, "Property name could not be found").ThrowAsJavaScriptException();
	return Napi::Value();
}

Napi::Int32Array Engine::CreateBufferedFacesArray(Napi::Env env, const Eigen::MatrixXi& F) const
{
	const uint32_t entries_per_face = 3;
	auto buffered_faces_array = Napi::Int32Array::New(env, entries_per_face * F.rows());

	//#pragma omp parallel for
	for (int32_t face_index = 0; face_index < F.rows(); face_index++)
	{
		const int base_index = entries_per_face * face_index;
		auto face = F.row(face_index);
		for (uint32_t i = 0; i < 3; i++)
		{
			buffered_faces_array[base_index + i] = face.coeffRef(i);
		}
	}

	return buffered_faces_array;
}

Napi::Int32Array Engine::CreateBufferedEdgesArray(Napi::Env env, const Eigen::MatrixXi& E) const
{
	const uint32_t entries_per_edge = 2;
	auto buffered_edges_array = Napi::Int32Array::New(env, entries_per_edge * E.rows());

	//#pragma omp parallel for
	for (int32_t edge_index = 0; edge_index < E.rows(); edge_index++)
	{
		const int base_index = entries_per_edge * edge_index;
		auto edge = E.row(edge_index);
		for (uint32_t i = 0; i < entries_per_edge; i++)
		{
			buffered_edges_array[base_index + i] = edge.coeffRef(i);
		}
	}

	return buffered_edges_array;
}

Napi::Value Engine::NativeToJS(Napi::Env env, const std::any& property_value)
{
	if (property_value.type() == typeid(double))
	{
		return NativeToJS(env, std::any_cast<const double&>(property_value));
	}

	if (property_value.type() == typeid(Eigen::VectorXd))
	{
		return NativeToJS(env, std::any_cast<const Eigen::VectorXd&>(property_value));
	}

	if (property_value.type() == typeid(std::vector<RDS::VertexIndex>))
	{
		return NativeToJS(env, std::any_cast<const std::vector<RDS::VertexIndex>&>(property_value));
	}

	if (property_value.type() == typeid(std::string))
	{
		return NativeToJS(env, std::any_cast<const std::string&>(property_value));
	}

	Napi::TypeError::New(env, "Type not supported").ThrowAsJavaScriptException();
	return Napi::Value();
}

Napi::Value Engine::NativeToJS(Napi::Env env, const Eigen::VectorXd& property_value)
{
	const auto rows_count = property_value.rows();
	Napi::Array array = Napi::Array::New(env, rows_count);
	for (int32_t i = 0; i < rows_count; i++)
	{
		array[i] = Napi::Number::New(env, property_value.coeffRef(i));
	}

	return array;
}

Napi::Value Engine::NativeToJS(Napi::Env env, const std::vector<RDS::VertexIndex>& property_value)
{
	const auto vector_size = property_value.size();
	Napi::Array array = Napi::Array::New(env, vector_size);
	for (int32_t i = 0; i < vector_size; i++)
	{
		array[i] = Napi::Number::New(env, property_value.at(i));
	}

	return array;
}

Napi::Value Engine::NativeToJS(Napi::Env env, const double property_value)
{
	return Napi::Number::New(env, property_value);
}

Napi::Value Engine::NativeToJS(Napi::Env env, const std::string& property_value)
{
	return Napi::String::New(env, property_value);
}

std::any Engine::JSToNative(Napi::Env env, const Napi::Value& value)
{
	if (value.IsString())
	{
		return std::make_any<std::string>(value.ToString());
	}

	if (value.IsNumber())
	{
		return std::make_any<double>(value.ToNumber());
	}

	if (value.IsArray())
	{
		auto array = value.As<Napi::Array>();
		std::vector<double> double_vector;
		for(int i = 0; i < array.Length(); i++)
		{
			Napi::Value value = array[i];	
			double_vector.push_back(value.ToNumber());
		}
		return std::make_any<std::vector<double>>(double_vector);
	}

	Napi::TypeError::New(env, "Type not supported").ThrowAsJavaScriptException();
	return std::any();
}

Napi::Array Engine::CreateEdges(Napi::Env env, const Eigen::MatrixX2i& E)
{
	Napi::Array edges_array = Napi::Array::New(env, E.rows());
	for (int edge_index = 0; edge_index < E.rows(); edge_index++)
	{
		Napi::Object edge_object = Napi::Object::New(env);
		int e0 = E(edge_index, 0);
		int e1 = E(edge_index, 1);

		edge_object.Set("e0", e0);
		edge_object.Set("e1", e1);
		edges_array[edge_index] = edge_object;
	}

	return edges_array;
}

Napi::Array Engine::CreateFaces(Napi::Env env, const Eigen::MatrixX3i& F)
{
	Napi::Array faces_array = Napi::Array::New(env, F.rows());
	for (int face_index = 0; face_index < F.rows(); face_index++)
	{
		Napi::Object face_object = Napi::Object::New(env);
		int v0 = F(face_index, 0);
		int v1 = F(face_index, 1);
		int v2 = F(face_index, 2);

		face_object.Set("v0", v0);
		face_object.Set("v1", v1);
		face_object.Set("v2", v2);
		faces_array[face_index] = face_object;
	}

	return faces_array;
}

void Engine::TryUpdateImageVertices()
{
	Eigen::VectorXd approximation_vector;
	if (newton_method_->GetApproximation(approximation_vector))
	{
		auto image_vertices = Eigen::Map<const Eigen::MatrixX2d>(approximation_vector.data(), approximation_vector.rows() >> 1, 2);
		mesh_wrapper_->SetImageVertices(image_vertices);
	}
}

void Engine::SetPositionWeight(const Napi::CallbackInfo& info, const Napi::Value& value)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (!value.IsNumber())
	{
		Napi::TypeError::New(env, "value is expected to be a Number").ThrowAsJavaScriptException();
		return;
	}

	/**
	 * Set position weight
	 */
	Napi::Number number = value.As<Napi::Number>();
	double position_weight = number.DoubleValue();
	if (position_) {
		position_->SetWeight(position_weight);
	}
}

Napi::Value Engine::GetPositionWeight(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number position_weight = Napi::Number::New(env, position_->GetWeight());

	return position_weight;
}

void Engine::SetSeamlessWeight(const Napi::CallbackInfo& info, const Napi::Value& value)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (!value.IsNumber())
	{
		Napi::TypeError::New(env, "value is expected to be a Number").ThrowAsJavaScriptException();
		return;
	}

	/**
	 * Set seamless weight
	 */
	Napi::Number number = value.As<Napi::Number>();
	double seamless_weight = number.DoubleValue();
}

Napi::Value Engine::GetSeamlessWeight(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	return env.Null();
}

void Engine::SetLambda(const Napi::CallbackInfo& info, const Napi::Value& value)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (!value.IsNumber())
	{
		Napi::TypeError::New(env, "value is expected to be a Number").ThrowAsJavaScriptException();
		return;
	}

	/**
	 * Set lambda
	 */
	Napi::Number number = value.As<Napi::Number>();
	double lambda = number.DoubleValue();

	separation_->SetWeight(lambda);
	symmetric_dirichlet_->SetWeight(1.0 - lambda);
}

Napi::Value Engine::GetLambda(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number lambda = Napi::Number::New(env, separation_->GetWeight());

	return lambda;
}

void Engine::SetDelta(const Napi::CallbackInfo& info, const Napi::Value& value)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (!value.IsNumber())
	{
		Napi::TypeError::New(env, "value is expected to be a Number").ThrowAsJavaScriptException();
		return;
	}

	/**
	 * Set delta
	 */
	Napi::Number number = value.As<Napi::Number>();
	double delta = number.DoubleValue();
	separation_->SetDelta(delta);
}

Napi::Value Engine::GetDelta(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	Napi::Number delta = Napi::Number::New(env, separation_->GetDelta());

	return delta;
}

Napi::Value Engine::CreateObjectiveFunctionDataObject(Napi::Env env, std::shared_ptr<ObjectiveFunction<Eigen::StorageOptions::RowMajor, Eigen::VectorXd>> objective_function) const
{
	Napi::Object data_object = Napi::Object::New(env);
	Napi::Object data_object_internal = Napi::Object::New(env);
	Napi::Array value_per_vertex_array = Napi::Array::New(env, mesh_wrapper_->GetImageVerticesCount());

	data_object.Set("name", objective_function->GetName());
	data_object.Set("data", data_object_internal);
	data_object_internal.Set("value", objective_function->GetValue());
	data_object_internal.Set("gradientNorm", objective_function->GetGradient().norm());

	auto value_per_vertex = objective_function->GetValuePerVertex();
	for (int32_t i = 0; i < value_per_vertex.rows(); i++)
	{
		value_per_vertex_array[i] = Napi::Number::New(env, value_per_vertex.coeffRef(i));
	}

	return data_object;
}

Napi::Value Engine::GetObjectiveFunctionsData(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	auto objective_functions_count = summation_objective_->GetObjectiveFunctionsCount();
	Napi::Array objective_functions_data_array = Napi::Array::New(env, objective_functions_count + 1);
	for (std::uint32_t index = 0; index < summation_objective_->GetObjectiveFunctionsCount(); index++)
	{
		objective_functions_data_array[index] = CreateObjectiveFunctionDataObject(env, summation_objective_->GetObjectiveFunction(index));
	}

	objective_functions_data_array[objective_functions_count] = CreateObjectiveFunctionDataObject(env, summation_objective_);

	return objective_functions_data_array;
}

Napi::Value Engine::LoadModel(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (info.Length() >= 1)
	{
		if (!info[0].IsString())
		{
			Napi::TypeError::New(env, "First argument is expected to be a String").ThrowAsJavaScriptException();
			return Napi::Value();
		}
	}
	else
	{
		Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Load model
	 */
	Napi::String value = info[0].As<Napi::String>();
	std::string model_file_path = std::string(value);
	mesh_wrapper_->LoadModel(model_file_path);

	return env.Null();
}

Engine::ModelFileType Engine::GetModelFileType(std::string modelFilePath)
{
	std::string fileExtension = modelFilePath.substr(modelFilePath.find_last_of(".") + 1);
	transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);
	if (fileExtension == "obj")
	{
		return Engine::ModelFileType::OBJ;
	}
	else if (fileExtension == "off")
	{
		return Engine::ModelFileType::OFF;
	}

	return Engine::ModelFileType::UNKNOWN;
}

Napi::Value Engine::ResumeSolver(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	if (newton_method_) 
	{
		newton_method_->Resume();
	}

	return env.Null();
}

Napi::Value Engine::PauseSolver(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	if (newton_method_) 
	{
		newton_method_->Pause();
	}

	return env.Null();
}

Engine::AlgorithmType Engine::StringToAlgorithmType(const std::string& algorithm_type_string)
{
	std::string mutable_string = algorithm_type_string;
	std::transform(mutable_string.begin(), mutable_string.end(), mutable_string.begin(), ::tolower);
	if(mutable_string == "autocuts")
	{
		return Engine::AlgorithmType::AUTOCUTS;
	}
	else if (mutable_string == "autoquads")
	{
		return Engine::AlgorithmType::AUTOQUADS;
	}

	throw std::exception("Unknown algorithm type");
}

Napi::Value Engine::SetAlgorithmType(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);
	
	if (newton_method_ != nullptr)
	{
		/**
		 * Validate input arguments
		 */
		if (info.Length() >= 1)
		{
			if (!info[0].IsString())
			{
				Napi::TypeError::New(env, "First argument is expected to be a String").ThrowAsJavaScriptException();
				return Napi::Value();
			}
		}
		else
		{
			Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
			return Napi::Value();
		}

		/**
		 * Set algorithm type
		 */
		AlgorithmType algorithm_type = StringToAlgorithmType(info[0].ToString());

		switch (algorithm_type)
		{
		case AlgorithmType::AUTOCUTS:
			summation_objective_ = autocuts_summation_objective_;
			break;
		case AlgorithmType::AUTOQUADS:
			summation_objective_ = autoquads_summation_objective_;
			break;
		}

		Eigen::VectorXd x0 = newton_method_->GetX();
		newton_method_->Terminate();
		newton_method_.release();
		newton_method_ = std::make_unique<NewtonMethod<PardisoSolver, Eigen::StorageOptions::RowMajor>>(summation_objective_, x0);
		newton_method_->EnableFlipAvoidingLineSearch(mesh_wrapper_->GetImageFaces());
		newton_method_->Start();
	}
	return env.Null();
}

Napi::Value Engine::ConstrainFacePosition(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (info.Length() >= 1)
	{
		if (!info[0].IsNumber())
		{
			Napi::TypeError::New(env, "First argument is expected to be a Number").ThrowAsJavaScriptException();
			return Napi::Value();
		}
	}
	else
	{
		Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Create a new barycenter position constraint for the given face
	 */
	int64_t face_index = info[0].As<Napi::Number>().Int64Value();
	auto face = mesh_wrapper_->GetImageFaceVerticesIndicesSTL(face_index);
	auto V_im = mesh_wrapper_->GetImageVertices();
	Eigen::Vector2d barycenter = Utils::CalculateBarycenter(face, V_im);
	auto face_data_provider = face_to_face_data_provider_map_.at(face);
	auto barycenter_position_objective = std::make_shared<FaceBarycenterPositionObjective<Eigen::StorageOptions::RowMajor>>(mesh_wrapper_, face_data_provider, barycenter);
	position_->AddObjectiveFunction(barycenter_position_objective);
	summation_objective_->Initialize();
	face_to_position_objective_map_.insert(std::make_pair(face, barycenter_position_objective));

	return Napi::Value();
}

Napi::Value Engine::UpdateConstrainedFacePosition(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (info.Length() >= 3)
	{
		if (!info[0].IsNumber())
		{
			Napi::TypeError::New(env, "First argument is expected to be a Number").ThrowAsJavaScriptException();
			return Napi::Value();
		}

		if (!info[1].IsNumber())
		{
			Napi::TypeError::New(env, "Second argument is expected to be a Number").ThrowAsJavaScriptException();
			return Napi::Value();
		}

		if (!info[2].IsNumber())
		{
			Napi::TypeError::New(env, "Third argument is expected to be a Number").ThrowAsJavaScriptException();
			return Napi::Value();
		}
	}
	else
	{
		Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
		return Napi::Value();
	}
 	/**
	 * Move the barycenter constraint by the given offset
	 */
	int64_t face_index = info[0].As<Napi::Number>().Int64Value();
	double offset_x = info[1].As<Napi::Number>().DoubleValue();
	double offset_y = info[2].As<Napi::Number>().DoubleValue();
	
	Eigen::Vector2d offset = Eigen::Vector2d(offset_x, offset_y);
	RDS::Face face = mesh_wrapper_->GetImageFaceVerticesIndicesSTL(face_index);

	face_to_position_objective_map_.at(face)->MoveFacePosition(offset);

	return Napi::Value();
}

Napi::Value Engine::UnconstrainFacePosition(const Napi::CallbackInfo& info)
{
	Napi::Env env = info.Env();
	Napi::HandleScope scope(env);

	/**
	 * Validate input arguments
	 */
	if (info.Length() >= 1)
	{
		if (!info[0].IsNumber())
		{
			Napi::TypeError::New(env, "First argument is expected to be a Number").ThrowAsJavaScriptException();
			return Napi::Value();
		}
	}
	else
	{
		Napi::TypeError::New(env, "Invalid number of arguments").ThrowAsJavaScriptException();
		return Napi::Value();
	}

	/**
	 * Remove the barycenter constraint of the given face
	 */
	int64_t face_index = info[0].As<Napi::Number>().Int64Value();
	RDS::Face face = mesh_wrapper_->GetImageFaceVerticesIndicesSTL(face_index);

	position_->RemoveObjectiveFunction(face_to_position_objective_map_.at(face));
	summation_objective_->Initialize();
	face_to_position_objective_map_.erase(face);

	return Napi::Value();
}