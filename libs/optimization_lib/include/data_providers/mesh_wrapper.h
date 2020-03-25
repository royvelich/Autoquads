#pragma once
#ifndef OPTIMIZATION_LIB_MESH_WRAPPER_H
#define OPTIMIZATION_LIB_MESH_WRAPPER_H

// STL Includes
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <algorithm>
#include <functional>
#include <string>

// Boost includes
#include <boost/signals2/signal.hpp>
#include <boost/functional/hash.hpp>

// Eigen Includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Optimization lib includes
#include "../core/core.h"
#include "./mesh_data_provider.h"

class MeshWrapper : public MeshDataProvider
{
public:

	/**
	 * Public type definitions
	 */
	enum class SoupInitType {
		RANDOM,
		ISOMETRIC
	};

	using ModelLoadedCallback = void();

	using EV2EVMap = std::vector<RDS::EdgePairDescriptor>;
	using VI2VIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;
	using VI2FIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;
	using EI2FIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;
	using FI2VIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;
	using FI2EIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;
	using FI2FIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;
	
	/**
	 * Constructors and destructor
	 */
	MeshWrapper();
	MeshWrapper(const Eigen::MatrixX3d& v, const Eigen::MatrixX3i& f);
	MeshWrapper(const std::string& modelFilename);
	virtual ~MeshWrapper();

	/**
	 * Setters
	 */
	void SetImageVertices(const Eigen::MatrixX2d& v_im);

	/**
	 * Getters
	 */
	RDS::Faces GetImageFacesSTL() const;
	const Eigen::MatrixX2d& GetImageVertices() const;
	const Eigen::MatrixX2i& GetImageEdges() const;
	const VI2VIsMap& GetImageNeighbours() const;
	const VI2VIsMap& GetDomainVerticesToImageVerticesMap() const;
	
	const VI2FIsMap& GetDomainVertexFaceAdjacency() const;
	const EI2FIsMap& GetDomainEdgeFaceAdjacency() const;
	const FI2VIsMap& GetDomainFaceVertexAdjacency() const;
	const FI2EIsMap& GetDomainFaceEdgeAdjacency() const;

	const VI2FIsMap& GetImageVertexFaceAdjacency() const;
	const EI2FIsMap& GetImageEdgeFaceAdjacency() const;
	const FI2VIsMap& GetImageFaceVertexAdjacency() const;
	const FI2EIsMap& GetImageFaceEdgeAdjacency() const;

	/**
	 * Overrides
	 */
	const Eigen::MatrixX3i& GetImageFaces() const override;
	const Eigen::MatrixX3i& GetDomainFaces() const override;
	const Eigen::MatrixX3d& GetDomainVertices() const override;
	const Eigen::MatrixX2i& GetDomainEdges() const override;
	const Eigen::MatrixX3d& GetD1() const override;
	const Eigen::MatrixX3d& GetD2() const override;
	const Eigen::SparseMatrix<double>& GetCorrespondingVertexPairsCoefficients() const override;
	const Eigen::VectorXd& GetCorrespondingVertexPairsEdgeLength() const override;
	int64_t GetImageVerticesCount() const override;
	int64_t GetDomainVerticesCount() const override;
	int64_t GetImageEdgesCount() const override;
	int64_t GetDomainEdgesCount() const override;
	const RDS::EdgePairDescriptors& GetEdgePairDescriptors() const;
	const RDS::FaceFans& GetFaceFans() const override;
	RDS::SparseVariableIndex GetXVariableIndex(RDS::VertexIndex vertex_index) const override;
	RDS::SparseVariableIndex GetYVariableIndex(RDS::VertexIndex vertex_index) const override;
	RDS::VertexIndex GetVertexIndex(RDS::SparseVariableIndex variable_index) const override;
	RDS::VertexIndex GetDomainVertexIndex(RDS::VertexIndex image_vertex_index) const override;
	int64_t GetVariablesCount() const override;
	RDS::EdgeIndex GetDomainEdgeIndex(RDS::EdgeDescriptor image_edge_descriptor) const override;
	RDS::EdgeIndex GetImageEdgeIndex(RDS::EdgeDescriptor image_edge_descriptor) const override;
	RDS::EdgeIndices GetDomainAdjacentEdgeIndicesByVertex(RDS::VertexIndex vertex_index) const override;
	RDS::EdgeIndices GetImageAdjacentEdgeIndicesByVertex(RDS::VertexIndex vertex_index) const override;
	
	/**
	 * Public methods
	 */
	Eigen::VectorXi GetImageFaceVerticesIndices(int64_t face_index);
	RDS::Face GetImageFaceVerticesIndicesSTL(int64_t face_index);
	Eigen::MatrixX2d GetImageVertices(const Eigen::VectorXi& vertex_indices);
	void LoadModel(const std::string& model_file_path);
	void RegisterModelLoadedCallback(const std::function<ModelLoadedCallback>& model_loaded_callback);

private:
	/**
	 * Private type definitions
	 */
	using EdgeDescriptor = std::pair<int64_t, int64_t>;
	using ED2EIMap = std::unordered_map<EdgeDescriptor, int64_t, RDS::UnorderedPairHash, RDS::UnorderedPairEquals>;
	using VI2VIMap = std::unordered_map<int64_t, int64_t>;
	using EI2EIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;
	using EI2EIMap = std::unordered_map<int64_t, int64_t>;
	using VI2EIsMap = std::unordered_map<int64_t, std::vector<int64_t>>;

	/**
	 * Private functions
	 */
	void Initialize();
	
	/**
	* Private enums
	*/
	enum class ModelFileType
	{
		OBJ,
		OFF,
		UNKNOWN
	};

	/**
	 * General use mesh methods
	 */
	void ComputeEdges(const Eigen::MatrixX3i& f, Eigen::MatrixX2i& e);
	void NormalizeVertices(Eigen::MatrixX3d& v);

	/**
	 * Discrete operators
	 */
	void ComputeSurfaceGradientPerFace(const Eigen::MatrixX3d& v, const Eigen::MatrixX3i& f, Eigen::MatrixX3d& d1, Eigen::MatrixX3d& d2);

	/**
	 * Triangle soup methods
	 */

	// Soup generation
	void GenerateSoupFaces(const Eigen::MatrixX3i& f_in, Eigen::MatrixX3i& f_out);
	void FixFlippedFaces(const Eigen::MatrixX3i& f_im, Eigen::MatrixX2d& v_im);
	void GenerateRandom2DSoup(const Eigen::MatrixX3i& f_in, Eigen::MatrixX3i& f_out, Eigen::MatrixX2d& v_out);
	void GenerateIsometric2DSoup(const Eigen::MatrixX3i& f_in, const Eigen::MatrixX3d& v_in, const ED2EIMap& ed_2_ei, const EI2FIsMap& ei_dom_2_fi , Eigen::MatrixX3i& f_out, Eigen::MatrixX2d& v_out);
	void CalculateAxes(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, Eigen::Vector3d& axis0, Eigen::Vector3d& axis1);
	void ProjectVertexToPlane(const Eigen::Vector3d& v0_in, const Eigen::Vector3d& v1_in, const Eigen::Vector3d& v2_in, const Eigen::Vector3d& v0_out, const Eigen::Vector3d& v1_out, Eigen::Vector3d& v2_out);
	void GetOrderedVertices(RDS::EdgeIndex edge_index, RDS::FaceIndex face_index, Eigen::Vector3i& vertex_indices);
	
	// Edge descriptor -> edge index map
	void ComputeEdgeDescriptorMap(const Eigen::MatrixX2i& e, ED2EIMap& ed_2_ei);

	// Domain edge index <-> image edge index maps
	void ComputeEdgeIndexMaps();

	// Domain vertex index <-> image vertex index maps
	void ComputeVertexIndexMaps();

	// Vertex index -> edge indices maps
	void ComputeVertexToEdgeIndexMaps();

	// Face/edge/vertex adjacency maps
	void ComputeAdjacencyMaps(
		const Eigen::MatrixX3i& f,
		const ED2EIMap& ed_2_ei,
		VI2FIsMap& vi_2_fi,
		VI2FIsMap& ei_2_fi,
		VI2FIsMap& fi_2_vi,
		VI2FIsMap& fi_2_ei,
		FI2FIsMap& fi_2_fi);

	// Image vertices corresponding pairs and image edges corresponding pairs
	void ComputeCorrespondingPairs();
	void ComputeCorrespondingVertexPairsCoefficients();
	void ComputeCorrespondingVertexPairsEdgeLength();

	// Compute vertex neighbours
	void ComputeVertexNeighbours();

	// Compute adjacent faces vertices
	void ComputeFaceFans();

	/**
	 * Private methods
	 */
	ModelFileType GetModelFileType(const std::string& modelFilePath);

	/**
	 * Fields
	 */

	// Domain matrices
	Eigen::MatrixX3d v_dom_;
	Eigen::MatrixX3i f_dom_;
	Eigen::MatrixX2i e_dom_;

	// Image matrices
	Eigen::MatrixX2d v_im_;
	Eigen::MatrixX3i f_im_;
	Eigen::MatrixX2i e_im_;

	// Discrete partial-derivatives matrices
	Eigen::MatrixX3d d1_;
	Eigen::MatrixX3d d2_;

	// Image corresponding pairs
	std::vector<std::pair<int64_t, int64_t>> cv_pairs_;
	std::vector<std::pair<int64_t, int64_t>> ce_pairs_;
	RDS::EdgePairDescriptors edge_pair_descriptors_;
	Eigen::SparseMatrix<double> cv_pairs_coefficients_;
	Eigen::VectorXd cv_pairs_edge_length_;

	// Image neighbour vertices
	VI2VIsMap v_im_2_neighbours;

	// Faces
	RDS::FaceFans face_fans_;

	// Face fans
	RDS::Faces faces_;

	// Maps
	ED2EIMap ed_im_2_ei_im_;
	ED2EIMap ed_dom_2_ei_dom_;
	VI2VIsMap v_dom_2_v_im_;
	VI2VIMap v_im_2_v_dom_;
	EI2EIsMap e_dom_2_e_im_;
	EI2EIMap e_im_2_e_dom_;
	VI2EIsMap v_im_2_e_im_;

	VI2FIsMap vi_im_2_fi_im_;
	EI2FIsMap ei_im_2_fi_im_;
	FI2VIsMap fi_im_2_vi_im_;
	FI2EIsMap fi_im_2_ei_im_;
	FI2FIsMap fi_im_2_fi_im_;
	
	VI2FIsMap vi_dom_2_fi_dom_;
	EI2FIsMap ei_dom_2_fi_dom_;
	FI2VIsMap fi_dom_2_vi_dom_;
	FI2EIsMap fi_dom_2_ei_dom_;
	FI2FIsMap fi_dom_2_fi_dom_;
	
	// Boost signals
	boost::signals2::signal<ModelLoadedCallback> model_loaded_signal_;
};

#endif