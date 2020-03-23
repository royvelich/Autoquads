#pragma once

#include "libs/optimization_lib/include/objective_functions/ObjectiveFunction.h"

//#include <FEMSimLib/SimMeshElement.h>
//#include <FEMSimLib/Node.h>
//#include <MathLib/MathLib.h>
//#include <MathLib/Matrix.h>


typedef Eigen::Matrix<double, 2, 3> Matrix2x3;

typedef Eigen::AngleAxisd AngleAxisd;
typedef Eigen::VectorXd dVector;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::MatrixXd MatrixNxM;
typedef Eigen::Matrix3d Matrix3x3;
typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::Triplet<double> MTriplet;

typedef Eigen::Matrix2d Matrix2x2;
typedef Eigen::Matrix4d Matrix4x4;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector4d Vector4d;

/*
	This class implements Constant Strain Triangles elements in 3D
	Mostly the same as CSTElement2D except for matrix sizes
*/
class CSTriangle3D : public ObjectiveFunction {
	
private:
	//material parameters...
	double shearModulus = 50, bulkModulus = 50;
	//relates area/volume to the mass of the element
	double massDensity = 1;
	//keep track of the rest shape area
	Eigen::VectorXd restShapeArea;

	//the collection of nodes that define the triangle element
	//Node* n[3];
	//parameters needed for gradient and hessian of the energy
	//Eigen::Vector3d dEdx[3];
	Eigen::MatrixXd CurrV;
	std::vector < Matrix3x3> ddEdxdx[3][3];
	//tmp matrices used to speed up computation of the deformation gradient, green strain, etc
	std::vector<Matrix2x2> dXInv, strain;
	std::vector<Eigen::Matrix<double, 3, 2>> FF, dEdF;

	
	//as a deformation measure, we need to compute the deformation gradient F. F maps deformed vectors dx to undeformed coords dX: dx = F*dX.
	//for linear basis functions, an easy way to compute it is by looking at the matrix that maps deformed traingle/tet edges to their underformed counterparts (F = dx * inv(dX)).
	//void computeDeformationGradient(Eigen::Matrix<double, 3, 2>& dxdX, int fi);

	//sets important properties of the rest shape using the set of points passed in as parameters
	void setRestShapeFromCurrentConfiguration();

	void setYoungsModulusAndPoissonsRatio(double E, double nu) {
		shearModulus = E / (2 * (1 + nu));
		bulkModulus = E / (3 * (1 - 2 * nu));
	}

	void setShearModulusAndBulkModulus(double G, double K) {
		shearModulus = G;
		bulkModulus = K;
	}
	Eigen::VectorXd getMass();
	void addEnergyHessianTo(const dVector& x);

	virtual void init_hessian();
public:
	CSTriangle3D();
	~CSTriangle3D();
	virtual void init();
	
	virtual void updateX(const Eigen::VectorXd& X);
	virtual double value(const bool update);
	virtual void gradient(Eigen::VectorXd& g, const bool update);
	virtual void hessian();
};