#pragma once
#ifndef OPTIMIZATION_LIB_CORE_H
#define OPTIMIZATION_LIB_CORE_H

// STL includes
#include <unordered_map>

// Boost includes
#include <boost/container_hash/hash.hpp>

// Eigen includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Optimization lib includes
#include "./utils.h"

namespace RDS
{
	/**
	 * Index types
	 */
	using VertexIndex = std::int64_t;
	using EdgeIndex = std::int64_t;
	using FaceIndex = std::int64_t;
	using EdgeIndices = std::vector<EdgeIndex>;
	using EdgeDescriptor = std::pair<VertexIndex, VertexIndex>;
	using EdgePairDescriptor = std::pair<EdgeDescriptor, EdgeDescriptor>;
	using EdgePairDescriptors = std::vector<EdgePairDescriptor>;
	using SparseVariableIndex = std::int64_t;
	using DenseVariableIndex = std::int64_t;
	using HessianTripletIndex = std::int64_t;
	using PartialDerivativeIndex = std::int64_t;
	using HessianEntry = std::pair<DenseVariableIndex, DenseVariableIndex>;
	using Face = std::vector<VertexIndex>;
	using Faces = std::vector<Face>;
	using FaceFanSlice = std::pair<VertexIndex, std::pair<VertexIndex, VertexIndex>>;
	using FaceFan = std::vector<FaceFanSlice>;
	using FaceFans = std::vector<FaceFan>;
	using ProjectionDescriptor = std::pair<RDS::VertexIndex, Eigen::Vector2d>;

	/**
	 * Enums
	 */
	enum class CoordinateType : int32_t
	{
		X,
		Y
	};
	
	/**
	 * Custom hash and equals function objects for unordered_map
	 * https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
	 * https://stackoverflow.com/questions/35985960/c-why-is-boosthash-combine-the-best-way-to-combine-hash-values/35991300#35991300
	 */
	struct OrderedPairHash {
		template <class T1, class T2>
		std::size_t operator () (const std::pair<T1, T2>& pair) const
		{
			const auto first = static_cast<int64_t>(pair.first);
			const auto second = static_cast<int64_t>(pair.second);

			std::size_t seed = 0;

			boost::hash_combine(seed, first);
			boost::hash_combine(seed, second);
			return seed;
		}
	};

	struct OrderedPairEquals {
		template <class T1, class T2>
		bool operator () (const std::pair<T1, T2>& pair1, const std::pair<T1, T2>& pair2) const
		{
			const auto pair1_first = static_cast<int64_t>(pair1.first);
			const auto pair1_second = static_cast<int64_t>(pair1.second);
			const auto pair2_first = static_cast<int64_t>(pair2.first);
			const auto pair2_second = static_cast<int64_t>(pair2.second);

			return (pair1_first == pair2_first) && (pair1_second == pair2_second);
		}
	};

	struct UnorderedPairHash {
		template <class T1, class T2>
		std::size_t operator () (const std::pair<T1, T2>& pair) const
		{
			const auto first = static_cast<int64_t>(pair.first);
			const auto second = static_cast<int64_t>(pair.second);

			const auto minmax_pair = std::minmax(first, second);
			std::size_t seed = 0;

			boost::hash_combine(seed, minmax_pair.first);
			boost::hash_combine(seed, minmax_pair.second);
			return seed;
		}
	};

	struct UnorderedPairEquals {
		template <class T1, class T2>
		bool operator () (const std::pair<T1, T2>& pair1, const std::pair<T1, T2>& pair2) const
		{
			const auto pair1_first = static_cast<int64_t>(pair1.first);
			const auto pair1_second = static_cast<int64_t>(pair1.second);
			const auto pair2_first = static_cast<int64_t>(pair2.first);
			const auto pair2_second = static_cast<int64_t>(pair2.second);

			const auto minmax_pair1 = std::minmax(pair1_first, pair1_second);
			const auto minmax_pair2 = std::minmax(pair2_first, pair2_second);
			return (minmax_pair1.first == minmax_pair2.first) && (minmax_pair1.second == minmax_pair2.second);
		}
	};

	struct VectorHash {
		template <class T>
		std::size_t operator () (const std::vector<T>& vector) const
		{
			std::size_t seed = 0;
			std::vector<T> sorted_vector = vector;
			std::sort(sorted_vector.begin(), sorted_vector.end());

			for (auto value : sorted_vector)
			{
				boost::hash_combine(seed, value);
			}

			return seed;
		}

		std::size_t operator () (const Eigen::VectorXi& vector) const
		{
			const auto vector_internal = std::vector<int64_t>(vector.data(), vector.data() + vector.rows());
			return this->operator()(vector_internal);
		}
	};

	struct VectorEquals {
		template <class T>
		bool operator () (const std::vector<T>& vector1, const std::vector<T>& vector2) const
		{
			if (vector1.size() != vector2.size())
			{
				return false;
			}

			std::vector<T> sorted_vector1 = vector1;
			std::vector<T> sorted_vector2 = vector2;
			std::sort(sorted_vector1.begin(), sorted_vector1.end());
			std::sort(sorted_vector2.begin(), sorted_vector2.end());

			return sorted_vector1 == sorted_vector2;
		}

		bool operator () (const Eigen::VectorXi& vector1, const Eigen::VectorXi& vector2) const
		{
			const auto vector1_internal = std::vector<int64_t>(vector1.data(), vector1.data() + vector1.rows());
			const auto vector2_internal = std::vector<int64_t>(vector2.data(), vector2.data() + vector2.rows());
			return this->operator()(vector1_internal, vector2_internal);
		}
	};

	// https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
	struct MatrixBaseHash
	{
		template <class Derived>
		std::size_t operator () (const Eigen::MatrixBase<Derived>& matrix_base) const
		{
			return Utils::GenerateHash(matrix_base);
		}
	};

	struct EdgePairDescriptorHash
	{
		std::size_t operator () (const RDS::EdgePairDescriptor& edge_pair_descriptor) const
		{
			const auto edge_descriptor1 = static_cast<RDS::EdgeDescriptor>(edge_pair_descriptor.first);
			const auto edge_descriptor2 = static_cast<RDS::EdgeDescriptor>(edge_pair_descriptor.second);

			RDS::EdgeDescriptor minmax_edge_descriptor1 = std::minmax(edge_descriptor1.first, edge_descriptor1.second);
			RDS::EdgeDescriptor minmax_edge_descriptor2 = std::minmax(edge_descriptor2.first, edge_descriptor2.second);

			if ((minmax_edge_descriptor1.first > minmax_edge_descriptor2.first) ||
				(minmax_edge_descriptor1.first == minmax_edge_descriptor2.first) && (minmax_edge_descriptor1.second > minmax_edge_descriptor2.second))
			{
				std::swap(minmax_edge_descriptor1, minmax_edge_descriptor2);
			}

			std::size_t seed = 0;
			boost::hash_combine(seed, minmax_edge_descriptor1.first);
			boost::hash_combine(seed, minmax_edge_descriptor1.second);
			boost::hash_combine(seed, minmax_edge_descriptor2.first);
			boost::hash_combine(seed, minmax_edge_descriptor2.second);
			return seed;
		}
	};

	struct EdgePairDescriptorEquals
	{
		bool operator () (const RDS::EdgePairDescriptor& edge_pair_descriptor1, const RDS::EdgePairDescriptor& edge_pair_descriptor2) const
		{
			const auto edge_descriptor1_1 = static_cast<RDS::EdgeDescriptor>(edge_pair_descriptor1.first);
			const auto edge_descriptor2_1 = static_cast<RDS::EdgeDescriptor>(edge_pair_descriptor1.second);

			const auto edge_descriptor1_2 = static_cast<RDS::EdgeDescriptor>(edge_pair_descriptor2.first);
			const auto edge_descriptor2_2 = static_cast<RDS::EdgeDescriptor>(edge_pair_descriptor2.second);

			RDS::EdgeDescriptor minmax_edge_descriptor1_1 = std::minmax(edge_descriptor1_1.first, edge_descriptor1_1.second);
			RDS::EdgeDescriptor minmax_edge_descriptor2_1 = std::minmax(edge_descriptor2_1.first, edge_descriptor2_1.second);

			RDS::EdgeDescriptor minmax_edge_descriptor1_2 = std::minmax(edge_descriptor1_2.first, edge_descriptor1_2.second);
			RDS::EdgeDescriptor minmax_edge_descriptor2_2 = std::minmax(edge_descriptor2_2.first, edge_descriptor2_2.second);

			if ((minmax_edge_descriptor1_1.first > minmax_edge_descriptor2_1.first) ||
				(minmax_edge_descriptor1_1.first == minmax_edge_descriptor2_1.first) && (minmax_edge_descriptor1_1.second > minmax_edge_descriptor2_1.second))
			{
				std::swap(minmax_edge_descriptor1_1, minmax_edge_descriptor2_1);
			}

			if ((minmax_edge_descriptor1_2.first > minmax_edge_descriptor2_2.first) ||
				(minmax_edge_descriptor1_2.first == minmax_edge_descriptor2_2.first) && (minmax_edge_descriptor1_2.second > minmax_edge_descriptor2_2.second))
			{
				std::swap(minmax_edge_descriptor1_2, minmax_edge_descriptor2_2);
			}

			return (minmax_edge_descriptor1_1 == minmax_edge_descriptor1_2) && (minmax_edge_descriptor2_1 == minmax_edge_descriptor2_2);
		}
	};

	struct HessianEntryHash
	{
		std::size_t operator () (const RDS::HessianEntry& hessian_entry) const
		{
			const auto first = static_cast<int64_t>(hessian_entry.first);
			const auto second = static_cast<int64_t>(hessian_entry.second);

			std::size_t seed = 0;

			boost::hash_combine(seed, first);
			boost::hash_combine(seed, second);
			return seed;
		}
	};

	struct HessianEntryEquals
	{
		bool operator () (const RDS::HessianEntry& hessian_entry1, const RDS::HessianEntry& hessian_entry2) const
		{
			const auto hessian_entry1_first = static_cast<int64_t>(hessian_entry1.first);
			const auto hessian_entry1_second = static_cast<int64_t>(hessian_entry1.second);
			const auto hessian_entry2_first = static_cast<int64_t>(hessian_entry2.first);
			const auto hessian_entry2_second = static_cast<int64_t>(hessian_entry2.second);

			return (hessian_entry1_first == hessian_entry2_first) && (hessian_entry1_second == hessian_entry2_second);
		}
	};

	/**
	 * Map types
	 */
	using SparseVariableIndexToDenseVariableIndexMap = std::unordered_map<RDS::SparseVariableIndex, RDS::DenseVariableIndex>;
	using SparseVariableIndexToVertexIndexMap = std::unordered_map<RDS::SparseVariableIndex, RDS::VertexIndex>;
	using DenseVariableIndexToSparseVariableIndexMap = std::unordered_map<RDS::DenseVariableIndex, RDS::SparseVariableIndex>;
	using HessianEntryToTripletIndexMap = std::unordered_map<RDS::HessianEntry, RDS::HessianTripletIndex, HessianEntryHash, HessianEntryEquals>;

}

#endif