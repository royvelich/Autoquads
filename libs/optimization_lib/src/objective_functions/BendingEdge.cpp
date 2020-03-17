﻿#include "objective_functions/BendingEdge.h"

#define SGN(x) (((x)<0)?(-1):(1))

BendingEdge::BendingEdge() {
	name = "Bending Edge energy";
	w = 1;
	std::cout << name << " constructor" << std::endl;
}

BendingEdge::~BendingEdge() {
	std::cout << name << " destructor" << std::endl;
}

void BendingEdge::init()
{
	if (V.size() == 0 || F.size() == 0)
		throw name + " must define members V,F before init()!";

	calculateHinges();

	restAngle.resize(num_hinges);
	restEdgeLength.resize(num_hinges);
	restConst.resize(num_hinges);
	restArea.resize(num_hinges);
	angle.resize(num_hinges);
	x0.resize(num_hinges, 3);
	x1.resize(num_hinges, 3);
	x2.resize(num_hinges, 3);
	x3.resize(num_hinges, 3);

	setRestShapeFromCurrentConfiguration();
	init_hessian();
}

void BendingEdge::updateX(const Eigen::VectorXd& X)
{
	assert(X.rows() == (3 * V.rows()));
	for (int hi = 0; hi < num_hinges; hi++) {
		//X = [x,x, ... ,x,y,y, ... ,y,z,z, ... ,z]
		x0.row(hi) = Eigen::Vector3d(
			X(x0_index(hi) + (0 * V.rows())),	//x-coordinate
			X(x0_index(hi) + (1 * V.rows())),	//Y-coordinate
			X(x0_index(hi) + (2 * V.rows()))	//Z-coordinate
		); 
		x1.row(hi) = Eigen::Vector3d(
			X(x1_index(hi) + (0 * V.rows())),	//x-coordinate
			X(x1_index(hi) + (1 * V.rows())),	//Y-coordinate
			X(x1_index(hi) + (2 * V.rows()))	//Z-coordinate
		);
		x2.row(hi) = Eigen::Vector3d(
			X(x2_index(hi) + (0 * V.rows())),	//x-coordinate
			X(x2_index(hi) + (1 * V.rows())),	//Y-coordinate
			X(x2_index(hi) + (2 * V.rows()))	//Z-coordinate
		);
		x3.row(hi) = Eigen::Vector3d(
			X(x3_index(hi) + (0 * V.rows())),	//x-coordinate
			X(x3_index(hi) + (1 * V.rows())),	//Y-coordinate
			X(x3_index(hi) + (2 * V.rows()))	//Z-coordinate
		);
	}
	getAngle();
}

void BendingEdge::calculateHinges() {
	std::vector<std::vector<std::vector<int>>> TT;
	std::vector<Eigen::Vector2d> hinges;
	igl::triangle_triangle_adjacency(F, TT);
	assert(TT.size() == F.rows());
	
	///////////////////////////////////////////////////////////
	//Part 1 - Find unique hinges
	for (int fi = 0; fi < TT.size(); fi++) {
		std::vector< std::vector<int>> CurrFace = TT[fi];
		assert(CurrFace.size() == 3 && "Each face should be a triangle (not square for example)!");
		for (std::vector<int> hinge : CurrFace) {
			if (hinge.size() == 1) {
				//add this "hinge"
				int FaceIndex1 = fi;
				int FaceIndex2 = hinge[0];

				if (FaceIndex2 < FaceIndex1) {
					//Skip
					//This hinge already exists!
					//Empty on purpose
				}
				else {
					hinges.push_back(Eigen::Vector2d(FaceIndex1, FaceIndex2));
				}
			}
			else if (hinge.size() == 0) {
				//Skip
				//This triangle has no another adjacent triangle on that edge
				//Empty on purpose
			}
			else {
				//We shouldn't get here!
				//The mesh is invalid
				assert("Each triangle should have only one adjacent triangle on each edge!");
			}

		}
	}

	///////////////////////////////////////////////////////////
	//Part 2 - Find x0,x1,x2,x3 indecis for each hinge
	num_hinges = hinges.size();
	x0_index.resize(num_hinges);
	x1_index.resize(num_hinges);
	x2_index.resize(num_hinges);
	x3_index.resize(num_hinges);
	
	for (int hi = 0; hi < num_hinges; hi++) {
		//first triangle vertices
		int v1 = F(hinges[hi](0), 0);
		int v2 = F(hinges[hi](0), 1);
		int v3 = F(hinges[hi](0), 2);
		//second triangle vertices
		int V1 = F(hinges[hi](1), 0);
		int V2 = F(hinges[hi](1), 1);
		int V3 = F(hinges[hi](1), 2);

		/*
		* Here we should find x0,x1,x2,x3
		* from the given two triangles (v1,v2,v3),(V1,V2,V3)
		*
		*	x0--x2
		*  / \ /
		* x3--x1
		*
		*/
		if (v1 != V1 && v1 != V2 && v1 != V3) {
			x2_index(hi) = v1;
			x0_index(hi) = v2;
			x1_index(hi) = v3;

			if (V1 != v1 && V1 != v2 && V1 != v3)
				x3_index(hi) = V1;
			else if (V2 != v1 && V2 != v2 && V2 != v3)
				x3_index(hi) = V2;
			else
				x3_index(hi) = V3;
		}
		else if (v2 != V1 && v2 != V2 && v2 != V3) {
			x2_index(hi) = v2;
			x0_index(hi) = v1;
			x1_index(hi) = v3;

			if (V1 != v1 && V1 != v2 && V1 != v3)
				x3_index(hi) = V1;
			else if (V2 != v1 && V2 != v2 && V2 != v3)
				x3_index(hi) = V2;
			else
				x3_index(hi) = V3;
		}
		else {
			x2_index(hi) = v3;
			x0_index(hi) = v1;
			x1_index(hi) = v2;

			if (V1 != v1 && V1 != v2 && V1 != v3)
				x3_index(hi) = V1;
			else if (V2 != v1 && V2 != v2 && V2 != v3)
				x3_index(hi) = V2;
			else
				x3_index(hi) = V3;
		}
	}
}

void BendingEdge::setRestShapeFromCurrentConfiguration() {
	for (int hi = 0; hi < num_hinges; hi++) {
		Eigen::Vector3d x0 = V.row(x0_index(hi));// = n[0]->getWorldPosition();
		Eigen::Vector3d x1 = V.row(x1_index(hi));// = n[1]->getWorldPosition();
		Eigen::Vector3d x2 = V.row(x2_index(hi));// = n[2]->getWorldPosition();
		Eigen::Vector3d x3 = V.row(x3_index(hi));// = n[3]->getWorldPosition();
		Eigen::Vector3d e0 = x1 - x0;
		Eigen::Vector3d n1 = e0.cross(x2 - x0);
		Eigen::Vector3d n2 = (x3 - x0).cross(e0);
		int sign = SGN(n1.cross(n2).dot(e0));
		double l_n1 = n1.norm();
		double l_n2 = n2.norm();

		//update rest variables
		restAngle(hi) = acos(n1.dot(n2) / (l_n1*l_n2))*sign;
		restEdgeLength(hi) = e0.norm();
		restArea(hi) = 0.5 * (l_n1 + l_n2);
	}
	restConst = 3 * restEdgeLength.cwiseProduct(restEdgeLength).cwiseProduct(restArea.cwiseInverse());
}

void BendingEdge::getAngle() {
	for (int hi = 0; hi < num_hinges; hi++) {
		Eigen::Vector3d e0 = x1.row(hi) - x0.row(hi);
		Eigen::Vector3d n1 = e0.cross(x2.row(hi) - x0.row(hi));
		Eigen::Vector3d n2 = (x3.row(hi) - x0.row(hi)).cross(e0);
		int sign = SGN(n1.cross(n2).dot(e0));
		double ratio = std::max(std::min(n1.dot(n2) / (n1.norm()*n2.norm()), 1.0), -1.0);
		angle(hi) = acos(ratio)*sign;
	}
}

double BendingEdge::value(const bool update)
{
	Eigen::VectorXd d_angle = angle - restAngle;
	/* alternative energy:
	k*3*restEdgeLength*restEdgeLength/restArea  * (phi(angle) - phi(restAngle))^2
	phi(alpha) = 2 * tan(alpha/2)	and	psi(alpha) = tan(alpha/pow(2,n))	s.t. phi = 2^n psi
	*/
	Eigen::VectorXd Energy = k * d_angle.cwiseProduct(d_angle).cwiseProduct(restConst);

	double value = Energy.sum();
	
	if (update) {
		Efi.setZero();
		energy_value = value;
	}
	return value;
}

void BendingEdge::gradient(Eigen::VectorXd& g, const bool update)
{
	g.conservativeResize(V.rows() * 3);
	g.setZero();

	Eigen::VectorXd d_angle = angle - restAngle;
	// constant factor depending only on rest shape
	Eigen::VectorXd K = 2 * k * restConst;
	// combine with factor depending on current angle
	Eigen::VectorXd zeta = K.cwiseProduct(d_angle);// zeta changes when using alternative energy, see Hessian

	for (int hi = 0; hi < num_hinges; hi++) {
		Eigen::Vector3d e0 = x1.row(hi) - x0.row(hi);
		Eigen::Vector3d e1 = x2.row(hi) - x0.row(hi);
		Eigen::Vector3d e2 = x3.row(hi) - x0.row(hi);
		Eigen::Vector3d e3 = x2.row(hi) - x1.row(hi);
		Eigen::Vector3d e4 = x3.row(hi) - x1.row(hi);
		Eigen::Vector3d n1 = e0.cross(e1);
		Eigen::Vector3d n2 = e2.cross(e0);
		double l_e0 = e0.norm(); e0 /= l_e0;
		double l_e1 = e1.norm(); e1 /= l_e1;
		double l_e2 = e2.norm(); e2 /= l_e2;
		double l_e3 = e3.norm(); e3 /= l_e3;
		double l_e4 = e4.norm(); e4 /= l_e4;
		double l_n1 = n1.norm(); n1 /= l_n1;
		double l_n2 = n2.norm(); n2 /= l_n2;
		double angle_1 = acos(e0.dot(e1));
		double angle_2 = acos(e2.dot(e0));
		double angle_3 = acos(e3.dot(-e0));
		double angle_4 = acos((-e0).dot(e4));
		double h_1 = l_n1 / l_e1;
		double h_2 = l_n2 / l_e2;
		double h_3 = l_n1 / l_e3;
		double h_4 = l_n2 / l_e4;
		double h_01 = l_n1 / l_e0;
		double h_02 = l_n2 / l_e0;


		Eigen::Vector3d component;
		// dE/dx0
		component = (n1*cos(angle_3) / h_3 + n2 * cos(angle_4) / h_4)*zeta(hi);
		for (int i = 0; i < 3; ++i)
			g[x0_index(hi) + (i*V.rows())] += component[i];

		// dE/dx1
		component = (n1*cos(angle_1) / h_1 + n2 * cos(angle_2) / h_2)*zeta(hi);
		for (int i = 0; i < 3; ++i)
			g[x1_index(hi) + (i*V.rows())] += component[i];

		// dE/dx2
		component = n1 * (-zeta(hi)) / h_01;
		for (int i = 0; i < 3; ++i)
			g[x2_index(hi) + (i*V.rows())] += component[i];

		// dE/dx3
		component = n2 * (-zeta(hi)) / h_02;
		for (int i = 0; i < 3; ++i)
			g[x3_index(hi) + (i*V.rows())] += component[i];
	}

	if (update)
		gradient_norm = g.norm();
}

void BendingEdge::hessian() {
	// constant factors
	int index = 0;
	Eigen::VectorXd d_angle = angle - restAngle;
	Eigen::VectorXd K = 2 * k *restConst;
	Eigen::VectorXd zeta = K.cwiseProduct(d_angle);// for alternative enrgy: zeta = K*(phi-rest_phi)*(1+psi*psi)
	Eigen::VectorXd xi = K;// for alternative enrgy: xi = K*(1+psi*psi)*(2*(psi-rest_psi)*psi + (1+psi*psi))

	for (int hi = 0; hi < num_hinges; hi++) {
		//start copied code from gradient
		Eigen::Vector3d e0 = x1.row(hi) - x0.row(hi);
		Eigen::Vector3d e1 = x2.row(hi) - x0.row(hi);
		Eigen::Vector3d e2 = x3.row(hi) - x0.row(hi);
		Eigen::Vector3d e3 = x2.row(hi) - x1.row(hi);
		Eigen::Vector3d e4 = x3.row(hi) - x1.row(hi);
		Eigen::Vector3d n1 = e0.cross(e1);
		Eigen::Vector3d n2 = e2.cross(e0);
		double l_e0 = e0.norm(); e0 /= l_e0;
		double l_e1 = e1.norm(); e1 /= l_e1;
		double l_e2 = e2.norm(); e2 /= l_e2;
		double l_e3 = e3.norm(); e3 /= l_e3;
		double l_e4 = e4.norm(); e4 /= l_e4;
		double l_n1 = n1.norm(); n1 /= l_n1;
		double l_n2 = n2.norm(); n2 /= l_n2;
		double angle_1 = acos(e0.dot(e1));
		double angle_2 = acos(e2.dot(e0));
		double angle_3 = acos(e3.dot(-e0));
		double angle_4 = acos((-e0).dot(e4));
		double h_1 = l_n1 / l_e1;
		double h_2 = l_n2 / l_e2;
		double h_3 = l_n1 / l_e3;
		double h_4 = l_n2 / l_e4;
		double h_01 = l_n1 / l_e0;
		double h_02 = l_n2 / l_e0;
		//end copied code from gradient

		//Gradient of angle
		Eigen::Vector3d grad_angle[4];
		grad_angle[0] = n1 * cos(angle_3) / h_3 + n2 * cos(angle_4) / h_4;
		grad_angle[1] = n1 * cos(angle_1) / h_1 + n2 * cos(angle_2) / h_2;
		grad_angle[2] = -n1 / h_01;
		grad_angle[3] = -n2 / h_02;

		// vectors m
		Eigen::Vector3d m1 = -e1.cross(n1);
		Eigen::Vector3d m2 = e2.cross(n2);
		Eigen::Vector3d m3 = e3.cross(n1);
		Eigen::Vector3d m4 = -e4.cross(n2);
		Eigen::Vector3d m01 = e0.cross(n1);
		Eigen::Vector3d m02 = -e0.cross(n2);

		//Hessian of angle
		Eigen::Matrix3d H_angle[4][4];
		//Eigen is making v1 * v2.T harder than it should be
		Eigen::RowVector3d n1_t = n1.transpose();
		Eigen::RowVector3d n2_t = n2.transpose();
		Eigen::RowVector3d m1_t = m1.transpose();
		Eigen::RowVector3d m2_t = m2.transpose();
		Eigen::RowVector3d m3_t = m3.transpose();
		Eigen::RowVector3d m4_t = m4.transpose();
		Eigen::RowVector3d m01_t = m01.transpose();
		Eigen::RowVector3d m02_t = m02.transpose();
		Eigen::Matrix3d B1 = n1 * m01_t / (l_e0*l_e0);
		Eigen::Matrix3d B2 = n2 * m02_t / (l_e0*l_e0);

		H_angle[0][0] = cos(angle_3) / (h_3*h_3) * (m3 * n1_t + n1 * m3_t) - B1
			+ cos(angle_4) / (h_4*h_4) * (m4 * n2_t + n2 * m4_t) - B2;
		H_angle[0][1] = cos(angle_3) / (h_1*h_3) * m1*n1_t + cos(angle_1) / (h_1*h_3)*n1*m3_t + B1
			+ cos(angle_4) / (h_2*h_4)*m2*n2_t + cos(angle_2) / (h_2*h_4)*n2*m4_t + B2;
		H_angle[0][2] = cos(angle_3) / (h_3*h_01)*m01*n1_t - n1 * m3_t / (h_01*h_3);
		H_angle[0][3] = cos(angle_4) / (h_4*h_02)*m02*n2_t - n2 * m4_t / (h_02*h_4);
		H_angle[1][1] = cos(angle_1) / (h_1*h_1)*(m1*n1_t + n1 * m1_t) - B1
			+ cos(angle_2) / (h_2*h_2)*(m2*n2_t + n2 * m2_t) - B2;
		H_angle[1][2] = cos(angle_1) / (h_1*h_01)*m01*n1_t - n1 * m1_t / (h_01*h_1);
		H_angle[1][3] = cos(angle_2) / (h_2*h_02)*m02*n2_t - n2 * m2_t / (h_02*h_2);
		H_angle[2][2] = -(n1*m01_t + m01 * n1_t) / (h_01*h_01);
		H_angle[2][3].setZero();
		H_angle[3][3] = -(n2*m02_t + m02 * n2_t) / (h_02*h_02);
		for (int i = 1; i < 4; ++i)
			for (int j = i - 1; j >= 0; --j)
				H_angle[i][j] = H_angle[j][i].transpose();

		//H(E) =  zeta * H(angle) + xi * grad(angle).transpose * grad(angle)
		Eigen::RowVector3d grad_angle_t[4];
		for (int i = 0; i < 4; ++i) grad_angle_t[i] = grad_angle[i].transpose();
		Eigen::Matrix3d H[4][4];
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				H[i][j] = zeta(hi) * H_angle[i][j] + xi(hi) * grad_angle[i] * grad_angle_t[j];

		//Finally, convert to triplets
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int ii = 0; ii < 3; ++ii) {
					int global_j;
					if (i == 0)
						global_j = x0_index(hi) + (ii*V.rows());
					else if (i == 1)
						global_j = x1_index(hi) + (ii*V.rows());
					else if (i == 2)
						global_j = x2_index(hi) + (ii*V.rows());
					else if (i == 3)
						global_j = x3_index(hi) + (ii*V.rows());

					for (int jj = 0; jj < 3; ++jj) {
						int global_i;

						if (j == 0)
							global_i = x0_index(hi) + (jj*V.rows());
						else if (j == 1)
							global_i = x1_index(hi) + (jj*V.rows());
						else if (j == 2)
							global_i = x2_index(hi) + (jj*V.rows());
						else if (j == 3)
							global_i = x3_index(hi) + (jj*V.rows());

						if (global_i <= global_j) {
							//hesEntries.push_back(Eigen::Triplet<double>(global_i, global_j, H[i][j](ii, jj)));
							SS[index++] = H[i][j](ii, jj);
						}	
					}
				}
	}
}

void BendingEdge::init_hessian()
{
	II.clear();
	JJ.clear();
	
	for (int hi = 0; hi < num_hinges; hi++) {
		//Finally, convert to triplets
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int ii = 0; ii < 3; ++ii) {
					int global_j;
					if (i == 0)
						global_j = x0_index(hi) + (ii*V.rows());
					else if (i == 1)
						global_j = x1_index(hi) + (ii*V.rows());
					else if (i == 2)
						global_j = x2_index(hi) + (ii*V.rows());
					else if (i == 3)
						global_j = x3_index(hi) + (ii*V.rows());

					for (int jj = 0; jj < 3; ++jj) {
						int global_i;

						if (j == 0)
							global_i = x0_index(hi) + (jj*V.rows());
						else if (j == 1)
							global_i = x1_index(hi) + (jj*V.rows());
						else if (j == 2)
							global_i = x2_index(hi) + (jj*V.rows());
						else if (j == 3)
							global_i = x3_index(hi) + (jj*V.rows());

						if (global_i <= global_j) {
							II.push_back(global_i);
							JJ.push_back(global_j);
						}
					}
				}
	}

	SS = std::vector<double>(II.size(), 0.);
}
