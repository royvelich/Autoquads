#pragma once
#include "libs/optimization_lib/include/objective_functions/ObjectiveFunction.h"

/*
	This class implements Constant Strain Triangles elements in 3D
	Mostly the same as CSTElement2D except for matrix sizes
*/
class MembraneConstraints : public ObjectiveFunction {
	
private:
	//material parameters...
	double shearModulus, bulkModulus;
	//keep track of the rest shape area
	Eigen::VectorXd restShapeArea;
	Eigen::MatrixXd CurrV;
	//tmp matrices used to speed up computation of the deformation gradient, green strain, etc
	std::vector<Eigen::Matrix2d> dXInv, strain;
	std::vector<Eigen::Matrix<double, 3, 2>> FF;

	//sets important properties of the rest shape using the set of points passed in as parameters
	void setRestShapeFromCurrentConfiguration();
	virtual void init_hessian();
public:
	MembraneConstraints();
	~MembraneConstraints();
	virtual void init();
	
	virtual void updateX(const Eigen::VectorXd& X);
	virtual double value(const bool update);
	virtual void gradient(Eigen::VectorXd& g, const bool update);
	virtual void hessian();
};