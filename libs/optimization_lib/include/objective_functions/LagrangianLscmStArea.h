#pragma once
#include "libs/optimization_lib/include/objective_functions/ConstrainedObjectiveFunction.h"

class LagrangianLscmStArea : public ConstrainedObjectiveFunction
{	
public:
	LagrangianLscmStArea();
	virtual double value(const bool update) override;
	virtual void gradient(Eigen::VectorXd& g, const bool update) override;
	virtual void hessian() override;
	virtual double AugmentedValue(const bool update) override;

	double objectiveValue(const bool update);
	double lagrangianValue(const bool update);
	Eigen::VectorXd constrainedValue(const bool update);

	void lagrangianGradient(Eigen::VectorXd& g, const bool update);
	void AuglagrangGradWRTX(Eigen::VectorXd& g, const bool update);
	void objectiveGradient(Eigen::VectorXd& g, const bool update);
};