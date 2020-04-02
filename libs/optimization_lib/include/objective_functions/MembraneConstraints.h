#pragma once
#include "libs/optimization_lib/include/objective_functions/ObjectiveFunction.h"

/*
	This class implements Constant Strain Triangles elements in 3D
	Mostly the same as CSTElement2D except for matrix sizes
*/


class MembraneConstraints : public ObjectiveFunction {
	
private:
	// Jacobian of the parameterization per face
	Eigen::VectorXd a;
	Eigen::VectorXd b;
	Eigen::VectorXd c;
	Eigen::VectorXd d;
	Eigen::VectorXd detJ;
	
	//dense mesh derivative matrices
	Eigen::Matrix3Xd D1d, D2d;

	//material type
	Utils::Material type;
	//material parameters...
	double shearModulus, bulkModulus;
	//keep track of the rest shape area
	Eigen::VectorXd restShapeArea;
	Eigen::MatrixX3d CurrV;
	//tmp matrices used to speed up computation of the deformation gradient, green strain, etc
	std::vector<Eigen::Matrix2d> dXInv, strain;
	std::vector<Eigen::Matrix<double, 3, 2>> F;

	//sets important properties of the rest shape using the set of points passed in as parameters
	void setRestShapeFromCurrentConfiguration();
	virtual void init_hessian();
public:
	MembraneConstraints(Utils::Material type);
	~MembraneConstraints();
	virtual void init();
	
	virtual void updateX(const Eigen::VectorXd& X);
	virtual double value(const bool update);
	virtual void gradient(Eigen::VectorXd& g, const bool update);
	virtual void hessian();
};