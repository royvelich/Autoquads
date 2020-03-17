#pragma once

#include "libs/optimization_lib/include/solvers/solver.h"

class GradientDescentSolver : public solver
{
public:
	GradientDescentSolver(const int solverID) : solver(solverID) {}
	virtual void step() override {
		objective->updateX(X);
		currentEnergy = objective->value(true);
		objective->gradient(g, true);
		p = -g;
	}
	virtual bool test_progress() override {
		return true;
	}
	virtual void internal_init() override {
		objective->updateX(X);
		g.resize(X.size());
	}
};