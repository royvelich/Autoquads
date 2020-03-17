#include "solvers/EigenSparseSolver.h"
#include <vector>
#include <iostream>

using namespace std;

template <typename vectorTypeI, typename vectorTypeS>
EigenSparseSolver<vectorTypeI, vectorTypeS>::EigenSparseSolver()
{
}

template <typename vectorTypeI, typename vectorTypeS>
EigenSparseSolver<vectorTypeI, vectorTypeS>::~EigenSparseSolver()
{
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenSparseSolver<vectorTypeI, vectorTypeS>::set_pattern(const vectorTypeI &II, const vectorTypeI &JJ, const vectorTypeS &SS)
{
	perpareMatrix(II, JJ, SS);
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenSparseSolver<vectorTypeI, vectorTypeS>::analyze_pattern()
{
	solver.analyzePattern(full_A/*UpperTriangular_A*/);
	//assert(solver.info() == Eigen::Success && "analyzePattern failed!");
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenSparseSolver<vectorTypeI, vectorTypeS>::factorize(const vectorTypeI &II, const vectorTypeI &JJ, const vectorTypeS &SS)
{
	perpareMatrix(II,JJ,SS);
	solver.factorize(full_A);
	assert(solver.info() == Eigen::Success && "factorization failed!");
}

template <typename vectorTypeI, typename vectorTypeS>
Eigen::VectorXd EigenSparseSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd &rhs)
{
	Eigen::VectorXd x;
	x = solver.solve(rhs);
	MSE = (full_A * x - rhs).cwiseAbs2().sum();
	return x;
}

template <typename vectorTypeI, typename vectorTypeS>
void EigenSparseSolver<vectorTypeI, vectorTypeS>::perpareMatrix(const vectorTypeI &II, const vectorTypeI &JJ, const vectorTypeS &SS)
{
	Eigen::SparseMatrix<double> UpperTriangular_A = Utils::BuildMatrix(II,JJ,SS);
	full_A = UpperTriangular_A.selfadjointView<Eigen::Upper>();

	double min_eig_value = full_A.toDense().eigenvalues().real().minCoeff();
	cout << "before: min_eig_value = " << min_eig_value << endl;
	if (min_eig_value < 0) {
		for (int i = 0; i < full_A.rows(); i++) {
			full_A.coeffRef(i, i) = full_A.coeff(i, i) + (-min_eig_value + 1e-4);
		}
	}
	cout << "after: full min_eig_value = " << full_A.toDense().eigenvalues().real().minCoeff() << endl;
	full_A.makeCompressed();
}

template class EigenSparseSolver<std::vector<int, std::allocator<int> >, std::vector<double, std::allocator<double> > >;
