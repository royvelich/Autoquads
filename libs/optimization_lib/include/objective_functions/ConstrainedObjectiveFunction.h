#pragma once

#include "libs/optimization_lib/include/objective_functions/TriangleMeshObjectiveFunction.h"

class ConstrainedObjectiveFunction : public TriangleMeshObjectiveFunction
{	
protected:
	Eigen::VectorXd X;
	
	virtual bool update_variables(const Eigen::VectorXd& X) override;
	virtual void init_hessian() override;
	void init_aug_hessian();
	void AddElementToAugHessian(std::vector<int> ind);
public:
	Eigen::VectorXd lambda;

	ConstrainedObjectiveFunction() {}
	virtual void init() override;
	virtual double AugmentedValue(const bool update) override;
	virtual void gradient(Eigen::VectorXd& g, const bool update) override;
	double lagrangianValue(const bool update);
	virtual double value(const bool update) override;
	
	virtual double objectiveValue(const bool update) = 0;
	virtual Eigen::VectorXd objectiveGradient(const bool update) = 0;
	virtual void objectiveHessian(std::vector<int>& I, std::vector<int>& J, std::vector<double>& S) = 0;
	virtual Eigen::VectorXd constrainedValue(const bool update) = 0;
	virtual void constrainedGradient(std::vector<int>& I, std::vector<int>& J, std::vector<double>& S) = 0;
	virtual void constrainedHessian(std::vector<std::vector<int>>& Is, std::vector < std::vector<int>>& Js, std::vector < std::vector<double>>& Ss) = 0;
	
	void lagrangianGradient(Eigen::VectorXd& g, const bool update);
	void AuglagrangGradWRTX(Eigen::VectorXd& g, const bool update);
	void hessian();
	void aughessian();

	//save values to show results
	double objective_value = 0;
	double constraint_value = 0;
	double objective_gradient_norm = 0;
	double constraint_gradient_norm = 0;
	float augmented_value_parameter = 9999.0f;
	// Hessian sparse reprensentation
	std::vector<int> II_aug, JJ_aug;
	std::vector<double> SS_aug;
};
