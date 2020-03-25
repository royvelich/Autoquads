#include "objective_functions/MembraneConstraints.h"

MembraneConstraints::MembraneConstraints() {
	type = Material::SYMMETRIC_DIRICHLET;
	if (type == Material::STVK) {
		name = "STVK";
	}
	else if (type == Material::SYMMETRIC_DIRICHLET) {
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
	CurrV = Eigen::Map<const Eigen::MatrixXd>(X.data(), X.rows() / 3, 3);
	F.clear();
	strain.clear();
	for (int fi = 0; fi < restShapeF.rows(); fi++) {
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

double MembraneConstraints::value(const bool update) {
	
	Eigen::VectorXd Energy(restShapeF.rows());
	for (int fi = 0; fi < restShapeF.rows(); fi++) {
		if (type == Material::STVK) {
			Energy(fi) = shearModulus * strain[fi].squaredNorm();
			Energy(fi) += (bulkModulus / 2) * pow(strain[fi].trace(), 2);
		}
		else if (type == Material::SYMMETRIC_DIRICHLET) {
			Energy(fi) = 0.5 * (1 + 1/ pow(strain[fi].determinant(),2)) * strain[fi].squaredNorm();
		}
	}

	double total_energy = restShapeArea.transpose() * Energy;

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
		Eigen::Matrix<double, 6, 9> dF_dX;
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
		
		if (type == Material::STVK) {
			dE_dstrain <<
				2 * shearModulus*strain[fi](0, 0) + bulkModulus * strain[fi].trace(),
				2 * shearModulus*strain[fi](0, 1),
				2 * shearModulus*strain[fi](1, 0),
				2 * shearModulus*strain[fi](1, 1) + bulkModulus * strain[fi].trace();
			dE_dstrain *= restShapeArea[fi];
		}
		else if (type == Material::SYMMETRIC_DIRICHLET) {
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
		

		Eigen::Matrix<double, 1, 9> dE_dX = dE_dstrain*dstrain_dF*dF_dX;

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
	for (int fi = 0; fi < restShapeF.rows(); fi++) {
		Eigen::Matrix<double, 6, 9> dF_dX;
		dF_dX <<
			-dXInv[fi](0, 0) - dXInv[fi](1, 0), dXInv[fi](0, 0), dXInv[fi](1, 0), 0, 0, 0, 0, 0, 0,
			-dXInv[fi](0, 1) - dXInv[fi](1, 1), dXInv[fi](0, 1), dXInv[fi](1, 1), 0, 0, 0, 0, 0, 0,
			0, 0, 0, -dXInv[fi](0, 0) - dXInv[fi](1, 0), dXInv[fi](0, 0), dXInv[fi](1, 0), 0, 0, 0,
			0, 0, 0, -dXInv[fi](0, 1) - dXInv[fi](1, 1), dXInv[fi](0, 1), dXInv[fi](1, 1), 0, 0, 0,
			0, 0, 0, 0, 0, 0, -dXInv[fi](0, 0) - dXInv[fi](1, 0), dXInv[fi](0, 0), dXInv[fi](1, 0),
			0, 0, 0, 0, 0, 0, -dXInv[fi](0, 1) - dXInv[fi](1, 1), dXInv[fi](0, 1), dXInv[fi](1, 1);

		Eigen::Matrix<double, 4, 6> dstrain_dF;
		dstrain_dF <<
			F[fi](0, 0), 0, F[fi](1, 0), 0, F[fi](2, 0), 0,
			0.5*F[fi](0, 1), 0.5*F[fi](0, 0), 0.5*F[fi](1, 1), 0.5*F[fi](1, 0), 0.5*F[fi](2, 1), 0.5*F[fi](2, 0),
			0.5*F[fi](0, 1), 0.5*F[fi](0, 0), 0.5*F[fi](1, 1), 0.5*F[fi](1, 0), 0.5*F[fi](2, 1), 0.5*F[fi](2, 0),
			0, F[fi](0, 1), 0, F[fi](1, 1), 0, F[fi](2, 1);

		Eigen::Matrix<double, 1, 4> dE_dstrain;
		Eigen::Matrix<double, 4, 4> dE_dstraindstrain;

		if (type == Material::STVK) {
			dE_dstrain <<
				2 * shearModulus*strain[fi](0, 0) + bulkModulus * strain[fi].trace(),
				2 * shearModulus*strain[fi](0, 1),
				2 * shearModulus*strain[fi](1, 0),
				2 * shearModulus*strain[fi](1, 1) + bulkModulus * strain[fi].trace();
			dE_dstrain *= restShapeArea[fi];

			dE_dstraindstrain <<
				2 * shearModulus + bulkModulus, 0, 0, bulkModulus,
				0, 2 * shearModulus, 0, 0,
				0, 0, 2 * shearModulus, 0,
				bulkModulus, 0, 0, 2 * shearModulus + bulkModulus;
			dE_dstraindstrain *= restShapeArea[fi];

		}
		else if (type == Material::SYMMETRIC_DIRICHLET) {
			double det = strain[fi].determinant();
			double a = strain[fi](0, 0);
			double b = strain[fi](0, 1);
			double c = strain[fi](1, 0);
			double d = strain[fi](1, 1);
			double Fnorm = strain[fi].squaredNorm();
			dE_dstrain <<
				a + a / pow(det, 2) - d * Fnorm / pow(det, 3),
				b + b / pow(det, 2) + c * Fnorm / pow(det, 3),
				c + c / pow(det, 2) + b * Fnorm / pow(det, 3),
				d + d / pow(det, 2) - a * Fnorm / pow(det, 3);
			dE_dstrain *= restShapeArea[fi];

			double aa = 1
				+ (1 / pow(det, 2))
				- ((4 * a*d) / pow(det, 3))
				+ ((3 * pow(d, 2)*Fnorm) / pow(det, 4));

			double bb = 1
				+ (1 / pow(det, 2))
				+ ((4 * b*c) / pow(det, 3))
				+ ((3 * pow(c, 2)*Fnorm) / pow(det, 4));

			double cc = 1
				+ (1 / pow(det, 2))
				+ ((4 * b*c) / pow(det, 3))
				+ ((3 * pow(b, 2)*Fnorm) / pow(det, 4));

			double dd = 1
				+ (1 / pow(det, 2))
				- ((4 * a*d) / pow(det, 3))
				+ ((3 * pow(a, 2)*Fnorm) / pow(det, 4));

			double ab = (-3 * c*d*Fnorm)
				+ (2 * (a*c - b*d)*det);
			ab /= pow(det, 4);

			double ac = (-3 * b*d*Fnorm)
				+ (2 * (a*b - c*d)*det);
			ac /= pow(det, 4);

			double ad = (3 * a*d*Fnorm)
				- ((2 * pow(a, 2) + 2 * pow(d, 2) + Fnorm)*det);
			ad /= pow(det, 4);

			double bc = (3 * b*c*Fnorm)
				+ ((2 * pow(b, 2) + 2 * pow(c, 2) + Fnorm)*det);
			bc /= pow(det, 4);

			double bd = (-3 * a*c*Fnorm)
				+ (2 * (c*d - a*b)*det);
			bd /= pow(det, 4);

			double cd = (-3 * a*b*Fnorm)
				+ (2 * (b*d - a*c)*det);
			cd /= pow(det, 4);

			dE_dstraindstrain <<
				aa, ab, ac, ad,
				ab, bb, bc, bd,
				ac, bc, cc, cd,
				ad, bd, cd, dd;
			dE_dstraindstrain *= restShapeArea[fi];
		}

		

		Eigen::Matrix<double, 6, 6> ds_dFdF___dE_ds; // ds_dFdF * dE_ds
		ds_dFdF___dE_ds <<
			dE_dstrain[0]							, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, 0										, 0										, 0										, 0										,
			0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, dE_dstrain[3]							, 0										, 0										, 0										, 0										,
			0										, 0										, dE_dstrain[0]							, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, 0										, 0										,
			0										, 0										, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, dE_dstrain[3]							, 0										, 0										,
			0										, 0										, 0										, 0										, dE_dstrain[0]							, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	,
			0										, 0										, 0										, 0										, 0.5*dE_dstrain[1] + 0.5*dE_dstrain[2]	, dE_dstrain[3];

		
		Eigen::Matrix<double, 9, 9> dE_dXdX =
			dF_dX.transpose() * dstrain_dF.transpose() * dE_dstraindstrain * dstrain_dF * dF_dX 
			+ dF_dX.transpose() * ds_dFdF___dE_ds * dF_dX;


		for (int v1 = 0; v1 < 3; v1++) {
			for (int v2 = 0; v2 < 3; v2++) {
				for (int xyz1 = 0; xyz1 < 3; xyz1++) {
					for (int xyz2 = 0; xyz2 < 3; xyz2++) {
						int global_i = restShapeF(fi, v1) + (xyz1*restShapeV.rows());
						int global_j = restShapeF(fi, v2) + (xyz2*restShapeV.rows());
						if (global_i <= global_j) {
							II.push_back(global_i);
							JJ.push_back(global_j);
							SS.push_back(dE_dXdX(3*xyz1 + v1, 3*xyz2 + v2));
						}
					}
				}
			}
		}

		/*Eigen::Matrix<double, 2, 3> FT = FF[fi].transpose();
		Eigen::Matrix<double, 2, 2> E = 0.5 * (FT * FF[fi]);
		E(0, 0) -= 0.5; E(1, 1) -= 0.5;
		Eigen::Matrix<double, 2, 2> I;
		I.setZero();
		Eigen::Matrix3d ddEdxdx[3][3];
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
			{
				Eigen::Matrix<double, 3, 2> dFdXij;
				dFdXij.setZero();
				if (i > 0)
					dFdXij(j, i - 1) = 1;
				else
					dFdXij(j, 0) = dFdXij(j, 1) = -1;
				dFdXij = dFdXij * dXInv[fi];
				Eigen::Matrix<double, 2, 2> dEdXij = 0.5 * (dFdXij.transpose() * FF[fi] + FT * dFdXij);
				Eigen::Matrix<double, 3, 2> dPdXij = dFdXij * (2.0 * shearModulus * E + bulkModulus * E.trace() * I);
				dPdXij += FF[fi] * (2.0 * shearModulus * dEdXij + bulkModulus * dEdXij.trace() * I);
				Eigen::Matrix<double, 3, 2> dHdXij = restShapeArea[fi] * dPdXij * dXInv[fi].transpose();

				for (int ii = 0; ii < 2; ++ii)
					for (int jj = 0; jj < 3; ++jj)
						ddEdxdx[ii + 1][i](jj, j) = dHdXij(jj, ii);
				ddEdxdx[0][i](0, j) = -dHdXij(0, 1) - dHdXij(0, 0);
				ddEdxdx[0][i](1, j) = -dHdXij(1, 1) - dHdXij(1, 0);
				ddEdxdx[0][i](2, j) = -dHdXij(2, 1) - dHdXij(2, 0);
			}
		
		II.clear();
		JJ.clear();
		SS.clear();
		for (int v1 = 0; v1 < 3; v1++) {
			for (int v2 = 0; v2 < 3; v2++) {
				for (int xyz1 = 0; xyz1 < 3; xyz1++) {
					for (int xyz2 = 0; xyz2 < 3; xyz2++) {
						int global_i = restShapeF(fi, v1) + (xyz1*restShapeV.rows());
						int global_j = restShapeF(fi, v2) + (xyz2*restShapeV.rows());
						if (global_i <= global_j) {
							II.push_back(global_i);
							JJ.push_back(global_j);
							SS.push_back(ddEdxdx[v1][v2](xyz1, xyz2));
						}
					}
				}
			}
		}*/
	}
	
}

