#pragma once

#include "libs/optimization_lib/include/solvers/solver.h"
#include "libs/optimization_lib/include/solvers/EigenSparseSolver.h"

class NewtonSolver : public solver
{
public:
	NewtonSolver(const int solverID): solver(solverID) {}

	virtual void step() override {
		objective->updateX(X);
		currentEnergy = objective->value(true);
		objective->gradient(g, true);
		objective->hessian();
		eigen_solver->factorize(objective->II, objective->JJ, objective->SS);
		Eigen::VectorXd rhs = -g;
		p = eigen_solver->solve(rhs);
	}

	virtual bool test_progress() override {
		return true;
	}

	virtual void internal_init() override {
		bool needs_init = eigen_solver == nullptr;
		if (needs_init) {
			eigen_solver = std::make_unique<EigenSparseSolver<std::vector<int>, std::vector<double>>>();
		}
		objective->updateX(X);
		g.resize(X.size());
		objective->hessian();
		if (needs_init) {
			eigen_solver->set_pattern(objective->II, objective->JJ, objective->SS);
			eigen_solver->analyze_pattern();
		}
	}
	
	Eigen::SparseMatrix<double> get_Hessian() {
		return eigen_solver->full_A;
	}

	double get_MSE() {
		return eigen_solver->MSE;
	}

	void SwitchPositiveDefiniteChecker() {
		eigen_solver->CheckPositiveDefinite = !(eigen_solver->CheckPositiveDefinite);
	}
	
private:
	std::unique_ptr<
		EigenSparseSolver<std::vector<int>, 
		std::vector<double>>
		> eigen_solver = nullptr;
};