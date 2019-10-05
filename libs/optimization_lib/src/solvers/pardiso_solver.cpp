// Optimization lib includes
#include <solvers/pardiso_solver.h>
#include <utils/utils.h>

// MKL includes
#include "mkl_pardiso.h"

PardisoSolver::PardisoSolver()
{
	mtype_ = -2;		/* Real symmetric matrix */
	nrhs_ = 1;			/* Number of right hand sides. */

	/* -------------------------------------*/
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------*/
	for (MKL_INT i = 0; i < 64; i++)
	{
		iparm_[i] = 0;
	}

	iparm_[0] = 1;		/* No solver default */
	iparm_[1] = 2;		/* Fill-in reordering from METIS */
	iparm_[3] = 0;		/* No iterative-direct algorithm */
	iparm_[4] = 0;		/* No user fill-in reducing permutation */
	iparm_[5] = 0;		/* Write solution into x */
	iparm_[7] = 1;		/* Max numbers of iterative refinement steps */
	iparm_[9] = 13;		/* Perturb the pivot elements with 1E-13 */
	iparm_[10] = 1;		/* Use nonsymmetric permutation and scaling MPS */
	iparm_[12] = 0;		/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm_[13] = 0;		/* Output: Number of perturbed pivots */
	iparm_[17] = -1;	/* Output: Number of nonzeros in the factor LU */
	iparm_[18] = -1;	/* Output: Mflops for LU factorization */
	iparm_[19] = 0;		/* Output: Numbers of CG Iterations */
	iparm_[34] = 1;		/* PARDISO use C-style indexing for ia and ja arrays */
	maxfct_ = 1;		/* Maximum number of numerical factorizations. */
	mnum_ = 1;			/* Which factorization to use. */
	msglvl_ = 0;		/* Print statistical information in file */
	error_ = 0;			/* Initialize error flag */

	/* ----------------------------------------------------------------*/
	/* .. Initialize the internal solver memory pointer. This is only  */
	/*   necessary for the FIRST call of the PARDISO solver.           */
	/* ----------------------------------------------------------------*/
	for (MKL_INT i = 0; i < 64; i++)
	{
		pt_[i] = 0;
	}
}

PardisoSolver::~PardisoSolver()
{
	/* --------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------*/

	/* Release internal memory. */
	phase_ = -1;
	PARDISO(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n_, &ddum_, ia_, ja_, &idum_, &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &error_);
}

void PardisoSolver::Solve(const Eigen::SparseMatrix<double, Eigen::StorageOptions::ColMajor>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{

}

void PardisoSolver::Solve(const Eigen::SparseMatrix<double, Eigen::StorageOptions::RowMajor>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
	auto A_copy = A;

	#pragma omp parallel for
	for (Eigen::DenseIndex i = 0; i < A.rows(); i++)
	{
		A_copy.coeffRef(i, i) = A_copy.coeffRef(i, i) + 0;
	}

	A_copy.makeCompressed();

	/* Matrix data. */
	n_ = A_copy.rows();
	ia_ = A_copy.outerIndexPtr();
	ja_ = A_copy.innerIndexPtr();
	a_ = A_copy.valuePtr();

	/* --------------------------------------------------------------------*/
	/* .. Reordering and Symbolic Factorization. This step also allocates  */
	/*    all memory that is necessary for the factorization.              */
	/* --------------------------------------------------------------------*/
	phase_ = 11;
	pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n_, a_, ia_, ja_, &idum_, &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &error_);

	/* ----------------------------*/
	/* .. Numerical factorization. */
	/* ----------------------------*/
	phase_ = 22;
	pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n_, a_, ia_, ja_, &idum_, &nrhs_, iparm_, &msglvl_, &ddum_, &ddum_, &error_);

	/* -----------------------------------------------*/
	/* .. Back substitution and iterative refinement. */
	/* -----------------------------------------------*/
	phase_ = 33;
	pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase_, &n_, a_, ia_, ja_, &idum_, &nrhs_, iparm_, &msglvl_, const_cast<double*>(b.data()), const_cast<double*>(x.data()), &error_);
}