#include "objective_functions/MembraneConstraints.h"

MembraneConstraints::MembraneConstraints(Utils::Material type) {
	this->type = type;
	if (type == Utils::STVK) {
		name = "STVK";
	}
	else if (type == Utils::SYMMETRIC_DIRICHLET) {
		name = "Symmetric Dirichlet";
	}
	w = 1;
	std::cout << name << " constructor" << std::endl;
}

MembraneConstraints::~MembraneConstraints() {
	std::cout << name << " destructor" << std::endl;
}

void MembraneConstraints::init()
{
	if (restShapeV.size() == 0 || restShapeF.size() == 0)
		throw name + " must define members V,F before init()!";
	
	assert(restShapeV.col(2).isZero() && "Warning: Rest shape is assumed to be in the plane (z coordinate must be zero in the beginning)");
	shearModulus = 0.3;
	bulkModulus = 1.5;
	setRestShapeFromCurrentConfiguration();
	init_hessian();
}

void MembraneConstraints::setRestShapeFromCurrentConfiguration() {

	a.resize(restShapeF.rows());
	b.resize(restShapeF.rows());
	c.resize(restShapeF.rows());
	d.resize(restShapeF.rows());
	detJ.resize(restShapeF.rows());
	
	Eigen::MatrixX3d D1cols, D2cols;
	Utils::computeSurfaceGradientPerFace(restShapeV, restShapeF, D1cols, D2cols);
	D1d = D1cols.transpose();
	D2d = D2cols.transpose();

	//edge vectors
	dXInv.clear();
	for (int fi = 0; fi < restShapeF.rows(); fi++) {
		Eigen::VectorXd V1 = restShapeV.row(restShapeF(fi, 1)) - restShapeV.row(restShapeF(fi, 0));
		Eigen::VectorXd V2 = restShapeV.row(restShapeF(fi, 2)) - restShapeV.row(restShapeF(fi, 0));
		//matrix that holds three edge vectors
		Eigen::Matrix2d dX; 
		dX << 
			V1[0], V2[0], 
			V1[1], V2[1];
		dXInv.push_back(dX.inverse()); //TODO .inverse() is baaad
	}
	
	//compute the area for each triangle
	igl::doublearea(restShapeV, restShapeF, restShapeArea);
	restShapeArea /= 2;
}

void MembraneConstraints::updateX(const Eigen::VectorXd& X)
{
	assert(X.rows() == (3 * restShapeV.rows()));
	CurrV = Eigen::Map<const Eigen::MatrixX3d>(X.data(), X.rows() / 3, 3);
	F.clear();
	strain.clear();
	Utils::LocalBasis(CurrV, restShapeF, B1, B2);

	for (int fi = 0; fi < restShapeF.rows(); fi++) {
		////////////////////////////////////////////////
		Eigen::Vector3d Xi, Yi;
		Xi <<
			CurrV.row(restShapeF(fi, 0)) * B1.row(fi).transpose(),
			CurrV.row(restShapeF(fi, 1)) * B1.row(fi).transpose(),
			CurrV.row(restShapeF(fi, 2)) * B1.row(fi).transpose();
		Yi << 
			CurrV.row(restShapeF(fi, 0)) * B2.row(fi).transpose(),
			CurrV.row(restShapeF(fi, 1)) * B2.row(fi).transpose(),
			CurrV.row(restShapeF(fi, 2)) * B2.row(fi).transpose();
		
		Eigen::Vector3d Dx = D1d.col(fi);
		Eigen::Vector3d Dy = D2d.col(fi);
		//prepare jacobian		
		a(fi) = Dx.transpose() * Xi;
		b(fi) = Dx.transpose() * Yi;
		c(fi) = Dy.transpose() * Xi;
		d(fi) = Dy.transpose() * Yi;
		detJ(fi) = a(fi) * d(fi) - b(fi) * c(fi);

		//std::cout << "a = " << a(fi) << std::endl;
		//std::cout << "b = " << b(fi) << std::endl;
		//std::cout << "c = " << c(fi) << std::endl;
		//std::cout << "d = " << d(fi) << std::endl;
		//std::cout << "detJ = " << detJ(fi) << std::endl;
		////////////////////////////////////////////////


		//edge vectors
		Eigen::VectorXd v1 = CurrV.row(restShapeF(fi, 1)) - CurrV.row(restShapeF(fi, 0));
		Eigen::VectorXd v2 = CurrV.row(restShapeF(fi, 2)) - CurrV.row(restShapeF(fi, 0));
		Eigen::Matrix<double, 3, 2> dx;
		dx <<
			v1(0), v2(0),
			v1(1), v2(1),
			v1(2), v2(2);
		
		F.push_back(dx * dXInv[fi]);

		//compute the Green Strain = 1/2 * (F'F-I)
		strain.push_back(F[fi].transpose() * F[fi]);
		strain[fi](0, 0) -= 1; strain[fi](1, 1) -= 1;
		strain[fi] *= 0.5;
	}
}

Eigen::Matrix<double, 1, 9> MembraneConstraints::da_dX(int fi) {
	Eigen::Vector3d Dx = D1d.col(fi);
	Eigen::Vector3d Dy = D2d.col(fi);
	
	Eigen::Matrix<double, 1, 3> V0 = CurrV.row(restShapeF(fi, 0));
	Eigen::Matrix<double, 1, 3> V1 = CurrV.row(restShapeF(fi, 1));
	Eigen::Matrix<double, 1, 3> V2 = CurrV.row(restShapeF(fi, 2));

	Eigen::Matrix<double, 1, 9> addx, addy, addz;
	addx << B1(fi, 0)	, 0			, 0			, B1(fi, 1)	, 0			, 0			, B1(fi, 2)	, 0			, 0;
	addy << 0			, B1(fi, 0)	, 0			, 0			, B1(fi, 1)	, 0			, 0			, B1(fi, 2)	, 0;
	addz << 0			, 0			, B1(fi, 0)	, 0			, 0			, B1(fi, 1)	, 0			, 0			, B1(fi, 2);
	
	Eigen::Matrix<double, 3, 9> XX, db1_dX = dB1_dX(fi);
	XX <<
		(V0 * db1_dX + addx),
		(V1 * db1_dX + addy),
		(V2 * db1_dX + addz);
	
	Eigen::Matrix<double, 1, 9> da_dX = Dx.transpose()*XX;
	Eigen::Matrix<double, 1, 9> dc_dX = Dy.transpose()*XX;
	return da_dX;
}

Eigen::Matrix<double, 9, 9> MembraneConstraints::ddB1_dXdX(int fi) {
	Eigen::Matrix<double, 9, 9> g_x,g_y,g_z;

	Eigen::Matrix<double, 3, 1> V0 = CurrV.row(restShapeF(fi, 0));
	Eigen::Matrix<double, 3, 1> V1 = CurrV.row(restShapeF(fi, 1));
	Eigen::Matrix<double, 3, 1> V2 = CurrV.row(restShapeF(fi, 2));
	double norm_V1V0 = (V1 - V0).norm();
	double a_b = V1[0] - V0[0]; // x1 - x0
	double c_d = V1[1] - V0[1]; // y1 - y0
	double e_f = V1[2] - V0[2]; // z1 - z0


	double xxx = (-3 * a_b*(pow(c_d, 2) + pow(e_f, 2))) / pow(norm_V1V0, 5);
	double yyy = (-3 * c_d*(pow(a_b, 2) + pow(e_f, 2))) / pow(norm_V1V0, 5);
	double zzz = (-3 * e_f*(pow(a_b, 2) + pow(c_d, 2))) / pow(norm_V1V0, 5);
	
	double xxy = (c_d*(pow(norm_V1V0, 2) - 3 * pow(a_b, 2))) / pow(norm_V1V0, 5);
	double xxz = (e_f*(pow(norm_V1V0, 2) - 3 * pow(a_b, 2))) / pow(norm_V1V0, 5);
	double xyy = (a_b*(pow(norm_V1V0, 2) - 3 * pow(c_d, 2))) / pow(norm_V1V0, 5);
	double xzz = (a_b*(pow(norm_V1V0, 2) - 3 * pow(e_f, 2))) / pow(norm_V1V0, 5);
	double yyz = (e_f*(pow(norm_V1V0, 2) - 3 * pow(c_d, 2))) / pow(norm_V1V0, 5);
	double yzz = (c_d*(pow(norm_V1V0, 2) - 3 * pow(e_f, 2))) / pow(norm_V1V0, 5);

	double xyz = (3*a_b*c_d*e_f) / pow(norm_V1V0, 5);

	g_x <<
		xxx	, -xxx	, 0, -xxy	, xxy	, 0, -xxz	, xxz	, 0,
		-xxx, xxx	, 0, xxy	, -xxy	, 0, xxz	, -xxz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0,
		-xxy, xxy	, 0, -xyy	, xyy	, 0, xyz	, -xyz	, 0,
		xxy	, -xxy	, 0, xyy	, -xyy	, 0, -xyz	, xyz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0,
		-xxz, xxz	, 0, xyz	, -xyz	, 0, -xzz	, xzz	, 0,
		xxz	, -xxz	, 0, -xyz	, xyz	, 0, xzz	, -xzz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0;
	
	g_y <<
		-xxy, xxy	, 0, -xyy	, xyy	, 0, xyz	, -xyz	, 0,
		xxy	, -xxy	, 0, xyy	, -xyy	, 0, -xyz	, xyz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0,
		-xyy, xyy	, 0, yyy	, -yyy	, 0, -yyz	, yyz	, 0,
		xyy	, -xyy	, 0, -yyy	, yyy	, 0, yyz	, -yyz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0,
		xyz	, -xyz	, 0, -yyz	, yyz	, 0, -yzz	, yzz	, 0,
		-xyz, xyz	, 0, yyz	, -yyz	, 0, yzz	, -yzz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0;
	
	g_z <<
		-xxz, xxz	, 0, xyz	, -xyz	, 0, -xzz	, xzz	, 0,
		xxz	, -xxz	, 0, -xyz	, xyz	, 0, xzz	, -xzz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0,
		xyz	, -xyz	, 0, -yyz	, yyz	, 0, -yzz	, yzz	, 0,
		-xyz, xyz	, 0, yyz	, -yyz	, 0, yzz	, -yzz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0,
		-xzz, xzz	, 0, -yzz	, yzz	, 0, zzz	, -zzz	, 0,
		xzz	, -xzz	, 0, yzz	, -yzz	, 0, -zzz	, zzz	, 0,
		0	, 0		, 0, 0		, 0		, 0, 0		, 0		, 0;

	return g_z;
}

Eigen::Matrix<double, 3, 9> MembraneConstraints::dB1_dX(int fi) {
	Eigen::Matrix<double, 3, 9> g;
	
	Eigen::Matrix<double, 3, 1> V0 = CurrV.row(restShapeF(fi, 0));
	Eigen::Matrix<double, 3, 1> V1 = CurrV.row(restShapeF(fi, 1));
	Eigen::Matrix<double, 3, 1> V2 = CurrV.row(restShapeF(fi, 2));
	double norm_V1V0 = (V1 - V0).norm();
	double a_b = V1[0] - V0[0]; // x1 - x0
	double c_d = V1[1] - V0[1]; // y1 - y0
	double e_f = V1[2] - V0[2]; // z1 - z0	

	// dB1.x/dx0
	double dB1x_dx0 = -(pow(c_d,2) + pow(e_f, 2)) / pow(norm_V1V0, 3);
	// dB1.y/dy0
	double dB1y_dy0 = -(pow(a_b, 2) + pow(e_f, 2)) / pow(norm_V1V0, 3);
	// dB1.z/dz0
	double dB1z_dz0 = -(pow(a_b, 2) + pow(c_d, 2)) / pow(norm_V1V0, 3);
	// dB1.x/dy0
	double dB1x_dy0 = (c_d*a_b) / pow(norm_V1V0, 3);
	// dB1.x/dz0
	double dB1x_dz0 = (e_f*a_b) / pow(norm_V1V0, 3);
	// dB1.y/dz0
	double dB1y_dz0 = (e_f*c_d) / pow(norm_V1V0, 3);
	g <<
		dB1x_dx0, -dB1x_dx0, 0, dB1x_dy0, -dB1x_dy0, 0, dB1x_dz0, -dB1x_dz0, 0,
		dB1x_dy0, -dB1x_dy0, 0, dB1y_dy0, -dB1y_dy0, 0, dB1y_dz0, -dB1y_dz0, 0,
		dB1x_dz0, -dB1x_dz0, 0, dB1y_dz0, -dB1y_dz0, 0, dB1z_dz0, -dB1z_dz0, 0;
	
	return g;
}


double MembraneConstraints::value(const bool update) {
	
	Eigen::VectorXd Energy(restShapeF.rows());
	/*for (int fi = 0; fi < restShapeF.rows(); fi++) {
		if (type == Utils::STVK) {
			Energy(fi) = shearModulus * strain[fi].squaredNorm();
			Energy(fi) += (bulkModulus / 2) * pow(strain[fi].trace(), 2);
		}
		else if (type == Utils::SYMMETRIC_DIRICHLET) {
			Energy(fi) = 0.5 * (1 + 1/ pow(strain[fi].determinant(),2)) * strain[fi].squaredNorm();
		}
	}
	double total_energy = restShapeArea.transpose() * Energy;
	*/
	
	//Energy = a;
	
	double total_energy = 0;
	for (int fi = 0; fi < restShapeF.rows(); fi++) {
		//total_energy += B1(fi, 0) + B1(fi, 1) + B1(fi, 2);
		//total_energy += CurrV(restShapeF(fi, 0), 0) * B1(fi, 0);

		//total_energy += CurrV.row(restShapeF(fi, 0)) * B1.row(fi).transpose();
		//total_energy += CurrV.row(restShapeF(fi, 1)) * B1.row(fi).transpose();
		//total_energy += CurrV.row(restShapeF(fi, 2)) * B1.row(fi).transpose();
		
		//total_energy += a(fi);
		total_energy += B1(fi,2);
	}
	

	if (update) {
		Efi = Energy;
		energy_value = total_energy;
	}
	return total_energy;
}

void MembraneConstraints::gradient(Eigen::VectorXd& g, const bool update)
{
	g.conservativeResize(restShapeV.rows() * 3);
	g.setZero();

	//compute the gradient of the energy using the chain rule: dE/dx = dE/dF * dF/dx. dE/dF is the first Piola-Kirchoff stress sensor, for which nice expressions exist.
	//compute the deformation gradient
	for (int fi = 0; fi < restShapeF.rows(); fi++) {
		/*Eigen::Matrix<double, 6, 9> dF_dX;
		dF_dX <<
			-dXInv[fi](0, 0) - dXInv[fi](1, 0)	, dXInv[fi](0, 0)	, dXInv[fi](1, 0)	, 0									, 0					, 0					, 0									, 0					, 0					,
			-dXInv[fi](0, 1) - dXInv[fi](1, 1)	, dXInv[fi](0, 1)	, dXInv[fi](1, 1)	, 0									, 0					, 0					, 0									, 0					, 0					,
			0									,0					,0					,-dXInv[fi](0, 0) - dXInv[fi](1, 0)	, dXInv[fi](0, 0)	, dXInv[fi](1, 0)	, 0									, 0					, 0					, 
			0									,0					,0					,-dXInv[fi](0, 1) - dXInv[fi](1, 1)	, dXInv[fi](0, 1)	, dXInv[fi](1, 1)	, 0									, 0					, 0					,
			0									,0					,0					,0									,0					,0					,-dXInv[fi](0, 0) - dXInv[fi](1, 0)	, dXInv[fi](0, 0)	, dXInv[fi](1, 0)	,
			0									,0					,0					,0									,0					,0					,-dXInv[fi](0, 1) - dXInv[fi](1, 1)	, dXInv[fi](0, 1)	, dXInv[fi](1, 1);
		
		Eigen::Matrix<double, 4, 6> dstrain_dF;
		dstrain_dF <<
			F[fi](0, 0)	, 0					, F[fi](1, 0)		, 0					, F[fi](2, 0)		, 0					,
			0.5*F[fi](0, 1), 0.5*F[fi](0, 0)	, 0.5*F[fi](1, 1)	, 0.5*F[fi](1, 0)	, 0.5*F[fi](2, 1)	, 0.5*F[fi](2, 0)	,
			0.5*F[fi](0, 1), 0.5*F[fi](0, 0)	, 0.5*F[fi](1, 1)	, 0.5*F[fi](1, 0)	, 0.5*F[fi](2, 1)	, 0.5*F[fi](2, 0)	,
			0				, F[fi](0, 1)		, 0					, F[fi](1, 1)		, 0					, F[fi](2, 1);
		
		Eigen::Matrix<double, 1, 4> dE_dstrain;
		
		if (type == Utils::STVK) {
			dE_dstrain <<
				2 * shearModulus*strain[fi](0, 0) + bulkModulus * strain[fi].trace(),
				2 * shearModulus*strain[fi](0, 1),
				2 * shearModulus*strain[fi](1, 0),
				2 * shearModulus*strain[fi](1, 1) + bulkModulus * strain[fi].trace();
			dE_dstrain *= restShapeArea[fi];
		}
		else if (type == Utils::SYMMETRIC_DIRICHLET) {
			double det = strain[fi].determinant();
			double a = strain[fi](0, 0);
			double b = strain[fi](0, 1);
			double c = strain[fi](1, 0);
			double d = strain[fi](1, 1);
			double Fnorm = strain[fi].squaredNorm();
			dE_dstrain <<
				a + a / pow(det,2) - d * Fnorm / pow(det, 3),
				b + b / pow(det, 2) + c * Fnorm / pow(det, 3),
				c + c / pow(det, 2) + b * Fnorm / pow(det, 3),
				d + d / pow(det, 2) - a * Fnorm / pow(det, 3);
			dE_dstrain *= restShapeArea[fi];
		}
		Eigen::Matrix<double, 1, 9> dE_dX = dE_dstrain*dstrain_dF*dF_dX;*/


		//Eigen::Matrix<double, 1, 9> dE_dX = da_dX(fi);
		Eigen::Matrix<double, 1, 9> dE_dX = dB1_dX(fi).row(2);
			
		
		for (int vi = 0; vi < 3; vi++)
			for (int xyz = 0; xyz < 3; xyz++)
				g[restShapeF(fi, vi) + (xyz*restShapeV.rows())] += dE_dX[xyz*3 + vi];


		//Eigen::Matrix2d stress = 2 * shearModulus * strain[fi];
		//stress(0, 0) += bulkModulus * (strain[fi].trace());
		//stress(1, 1) += bulkModulus * (strain[fi].trace());

		////Eigen::Matrix<double, 3, 2> dEdF = FF[fi] * (stress + stress.transpose());
		//Eigen::Matrix<double, 3, 2> dEdF = FF[fi] * stress;

		////dF/dx is going to be some +/- Xinv terms. The forces on nodes 1,2 can be writen as: dE/dF * XInv', while the force on node 0 is -f1-f2;
		//Eigen::Vector3d dEdx[3];
		//dEdx[1] =
		//	Eigen::Vector3d(
		//		dEdF(0, 0) * dXInv[fi](0, 0) + dEdF(0, 1) * dXInv[fi](0, 1),
		//		dEdF(1, 0) * dXInv[fi](0, 0) + dEdF(1, 1) * dXInv[fi](0, 1),
		//		dEdF(2, 0) * dXInv[fi](0, 0) + dEdF(2, 1) * dXInv[fi](0, 1)
		//	) * restShapeArea[fi];

		//dEdx[2] =
		//	Eigen::Vector3d(
		//		dEdF(0, 0) * dXInv[fi](1, 0) + dEdF(0, 1) * dXInv[fi](1, 1),
		//		dEdF(1, 0) * dXInv[fi](1, 0) + dEdF(1, 1) * dXInv[fi](1, 1),
		//		dEdF(2, 0) * dXInv[fi](1, 0) + dEdF(2, 1) * dXInv[fi](1, 1)
		//	) * restShapeArea[fi];
		//dEdx[0] = -dEdx[1] - dEdx[2];

		//for (int vi = 0; vi < 3; vi++)
		//	for (int xyz = 0; xyz < 3; xyz++)
		//		g[restShapeF(fi, vi) + (xyz*restShapeV.rows())] += dEdx[vi](xyz);
	}

	if (update)
		gradient_norm = g.norm();
}

void MembraneConstraints::init_hessian() {

}

void MembraneConstraints::hessian() {
	II.clear();
	JJ.clear();
	SS.clear();
	for (int fi = 0; fi < restShapeF.rows(); fi++) {
		//Eigen::Matrix<double, 6, 9> dF_dX;
		//dF_dX <<
		//	-dXInv[fi](0, 0) - dXInv[fi](1, 0), dXInv[fi](0, 0), dXInv[fi](1, 0), 0, 0, 0, 0, 0, 0,
		//	-dXInv[fi](0, 1) - dXInv[fi](1, 1), dXInv[fi](0, 1), dXInv[fi](1, 1), 0, 0, 0, 0, 0, 0,
		//	0, 0, 0, -dXInv[fi](0, 0) - dXInv[fi](1, 0), dXInv[fi](0, 0), dXInv[fi](1, 0), 0, 0, 0,
		//	0, 0, 0, -dXInv[fi](0, 1) - dXInv[fi](1, 1), dXInv[fi](0, 1), dXInv[fi](1, 1), 0, 0, 0,
		//	0, 0, 0, 0, 0, 0, -dXInv[fi](0, 0) - dXInv[fi](1, 0), dXInv[fi](0, 0), dXInv[fi](1, 0),
		//	0, 0, 0, 0, 0, 0, -dXInv[fi](0, 1) - dXInv[fi](1, 1), dXInv[fi](0, 1), dXInv[fi](1, 1);

		//Eigen::Matrix<double, 4, 6> dstrain_dF;
		//dstrain_dF <<
		//	F[fi](0, 0), 0, F[fi](1, 0), 0, F[fi](2, 0), 0,
		//	0.5*F[fi](0, 1), 0.5*F[fi](0, 0), 0.5*F[fi](1, 1), 0.5*F[fi](1, 0), 0.5*F[fi](2, 1), 0.5*F[fi](2, 0),
		//	0.5*F[fi](0, 1), 0.5*F[fi](0, 0), 0.5*F[fi](1, 1), 0.5*F[fi](1, 0), 0.5*F[fi](2, 1), 0.5*F[fi](2, 0),
		//	0, F[fi](0, 1), 0, F[fi](1, 1), 0, F[fi](2, 1);

		//Eigen::Matrix<double, 1, 4> dE_dstrain;
		//Eigen::Matrix<double, 4, 4> dE_dstraindstrain;

		//if (type == Utils::STVK) {
		//	dE_dstrain <<
		//		2 * shearModulus*strain[fi](0, 0) + bulkModulus * strain[fi].trace(),
		//		2 * shearModulus*strain[fi](0, 1),
		//		2 * shearModulus*strain[fi](1, 0),
		//		2 * shearModulus*strain[fi](1, 1) + bulkModulus * strain[fi].trace();
		//	dE_dstrain *= restShapeArea[fi];

		//	dE_dstraindstrain <<
		//		2 * shearModulus + bulkModulus, 0, 0, bulkModulus,
		//		0, 2 * shearModulus, 0, 0,
		//		0, 0, 2 * shearModulus, 0,
		//		bulkModulus, 0, 0, 2 * shearModulus + bulkModulus;
		//	dE_dstraindstrain *= restShapeArea[fi];

		//}
		//else if (type == Utils::SYMMETRIC_DIRICHLET) {
		//	double det = strain[fi].determinant();
		//	double a = strain[fi](0, 0);
		//	double b = strain[fi](0, 1);
		//	double c = strain[fi](1, 0);
		//	double d = strain[fi](1, 1);
		//	double Fnorm = strain[fi].squaredNorm();
		//	dE_dstrain <<
		//		a + a / pow(det, 2) - d * Fnorm / pow(det, 3),
		//		b + b / pow(det, 2) + c * Fnorm / pow(det, 3),
		//		c + c / pow(det, 2) + b * Fnorm / pow(det, 3),
		//		d + d / pow(det, 2) - a * Fnorm / pow(det, 3);
		//	dE_dstrain *= restShapeArea[fi];

		//	double aa = 1
		//		+ (1 / pow(det, 2))
		//		- ((4 * a*d) / pow(det, 3))
		//		+ ((3 * pow(d, 2)*Fnorm) / pow(det, 4));

		//	double bb = 1
		//		+ (1 / pow(det, 2))
		//		+ ((4 * b*c) / pow(det, 3))
		//		+ ((3 * pow(c, 2)*Fnorm) / pow(det, 4));

		//	double cc = 1
		//		+ (1 / pow(det, 2))
		//		+ ((4 * b*c) / pow(det, 3))
		//		+ ((3 * pow(b, 2)*Fnorm) / pow(det, 4));

		//	double dd = 1
		//		+ (1 / pow(det, 2))
		//		- ((4 * a*d) / pow(det, 3))
		//		+ ((3 * pow(a, 2)*Fnorm) / pow(det, 4));

		//	double ab = (-3 * c*d*Fnorm)
		//		+ (2 * (a*c - b*d)*det);
		//	ab /= pow(det, 4);

		//	double ac = (-3 * b*d*Fnorm)
		//		+ (2 * (a*b - c*d)*det);
		//	ac /= pow(det, 4);

		//	double ad = (3 * a*d*Fnorm)
		//		- ((2 * pow(a, 2) + 2 * pow(d, 2) + Fnorm)*det);
		//	ad /= pow(det, 4);

		//	double bc = (3 * b*c*Fnorm)
		//		+ ((2 * pow(b, 2) + 2 * pow(c, 2) + Fnorm)*det);
		//	bc /= pow(det, 4);

		//	double bd = (-3 * a*c*Fnorm)
		//		+ (2 * (c*d - a*b)*det);
		//	bd /= pow(det, 4);

		//	double cd = (-3 * a*b*Fnorm)
		//		+ (2 * (b*d - a*c)*det);
		//	cd /= pow(det, 4);

		//	dE_dstraindstrain <<
		//		aa, ab, ac, ad,
		//		ab, bb, bc, bd,
		//		ac, bc, cc, cd,
		//		ad, bd, cd, dd;
		//	dE_dstraindstrain *= restShapeArea[fi];
		//}

		//

		//Eigen::Matrix<double, 6, 6> ds_dFdF___dE_ds; // ds_dFdF * dE_ds
		//ds_dFdF___dE_ds <<
		//	dE_dstrain[0]							, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, 0										, 0										, 0										, 0										,
		//	0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, dE_dstrain[3]							, 0										, 0										, 0										, 0										,
		//	0										, 0										, dE_dstrain[0]							, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, 0										, 0										,
		//	0										, 0										, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, dE_dstrain[3]							, 0										, 0										,
		//	0										, 0										, 0										, 0										, dE_dstrain[0]							, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	,
		//	0										, 0										, 0										, 0										, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, dE_dstrain[3];

		//
		//Eigen::Matrix<double, 9, 9> dE_dXdX =
		//	dF_dX.transpose() * dstrain_dF.transpose() * dE_dstraindstrain * dstrain_dF * dF_dX 
		//	+ dF_dX.transpose() * ds_dFdF___dE_ds * dF_dX;

		Eigen::Matrix<double, 9, 9> dE_dXdX = ddB1_dXdX(fi);

		for (int v1 = 0; v1 < 3; v1++) {
			for (int v2 = 0; v2 < 3; v2++) {
				for (int xyz1 = 0; xyz1 < 3; xyz1++) {
					for (int xyz2 = 0; xyz2 < 3; xyz2++) {
						int global_i = restShapeF(fi, v1) + (xyz1*restShapeV.rows());
						int global_j = restShapeF(fi, v2) + (xyz2*restShapeV.rows());
						//if (global_i <= global_j) {
							II.push_back(global_i);
							JJ.push_back(global_j);
							SS.push_back(dE_dXdX(3*xyz1 + v1, 3*xyz2 + v2));
						//}
					}
				}
			}
		}
	}
}

