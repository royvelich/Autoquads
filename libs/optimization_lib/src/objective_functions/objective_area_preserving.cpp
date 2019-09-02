#include <objective_functions/objective_area_preserving.h>
#include <limits>
#include <igl/doublearea.h>
//#include <chrono>

using namespace std;
using namespace Eigen;

ObjectiveAreaPreserving::ObjectiveAreaPreserving()
{
	name = "Objective: Area Preserving";
}
void ObjectiveAreaPreserving::init()
{
	if (V.size() == 0 || F.size() == 0)
		throw "DistortionAreaPreserving must define members V,F before init()!";

	numF = F.rows();
	numV = V.rows();

	Fuv.resize(6, numF);
	Fuv.topRows(3) = F.transpose();
	Fuv.bottomRows(3) = Fuv.topRows(3) + Eigen::MatrixXi::Constant(3, numF, numV);


	a.resize(numF);
	b.resize(numF);
	c.resize(numF);
	d.resize(numF);
	alpha.resize(numF, 2);
	beta.resize(numF, 2);
	s.resize(numF, 2);
	v.resize(numF, 4);
	u.resize(numF, 4);

	Dsd[0].resize(6, numF);
	Dsd[1].resize(6, numF);

	//Parameterization J mats resize
	detJ.resize(numF);
	grad.resize(numF, 6);

	// compute init energy matrices
	igl::doublearea(V, F, Area);
	Area /= 2;

	Eigen::MatrixX3d D1cols, D2cols;

	Eigen::MatrixX3d V3d;
	V3d.resize(V.rows(), 3);
	V3d.leftCols(2) = V;
	V3d.col(2).setZero();
	Utils::computeSurfaceGradientPerFace(V3d, F, D1cols, D2cols);
	D1d = D1cols.transpose();
	D2d = D2cols.transpose();

	//columns belong to different faces
	a1d.resize(6, numF);
	a2d.resize(6, numF);
	b1d.resize(6, numF);
	b2d.resize(6, numF);

	a1d.topRows(3) = 0.5*D1d;
	a1d.bottomRows(3) = 0.5*D2d;

	a2d.topRows(3) = -0.5*D2d;
	a2d.bottomRows(3) = 0.5*D1d;

	b1d.topRows(3) = 0.5*D1d;
	b1d.bottomRows(3) = -0.5*D2d;

	b2d.topRows(3) = 0.5*D2d;
	b2d.bottomRows(3) = 0.5*D1d;

	Hi.resize(numF);
	prepare_hessian();
	w = 1;
}

void ObjectiveAreaPreserving::updateX(const VectorXd& X)
{
	bool inversions_exist = updateJ(X);
}

double ObjectiveAreaPreserving::value()
{
	// E = (det(J) - 1)^2
	MatrixXd ones(detJ.rows(), 1);
	ones.setConstant(1);
	Efi = (detJ - ones).cwiseAbs2();

	double f = 0.5*(Area.asDiagonal()*Efi).sum();
	return f;
}

void ObjectiveAreaPreserving::gradient(VectorXd& g)
{
	UpdateSSVDFunction();
	g.conservativeResize(numV * 2);
	g.setZero();
	ComputeDenseSSVDDerivatives();

	for (int fi = 0; fi < numF; ++fi) {
		VectorXd gi;
		gi.resize(6);
		gi = Area(fi)*(grad.row(fi));
		for (int vi = 0; vi < 6; ++vi) {
			g(Fuv(vi, fi)) += gi(vi);
		}
	}
}

void ObjectiveAreaPreserving::hessian()
{
	UpdateSSVDFunction();
	ComputeDenseSSVDDerivatives();

	auto lambda1 = [](double a) {return a - 1.0 / (a*a*a); };
	//gradient of outer function in composition
	Eigen::VectorXd gradfS = s.col(0).unaryExpr(lambda1);
	Eigen::VectorXd gradfs = s.col(1).unaryExpr(lambda1);
	auto lambda2 = [](double a) {return 1 + 3 / (a*a*a*a); };
	//hessian of outer function in composition (diagonal)
	Eigen::VectorXd HS = s.col(0).unaryExpr(lambda2);
	Eigen::VectorXd Hs = s.col(1).unaryExpr(lambda2);
	//simliarity alpha
	Eigen::VectorXd aY = 0.5*(a + d);
	Eigen::VectorXd bY = 0.5*(c - b);
	//anti similarity beta
	Eigen::VectorXd cY = 0.5*(a - d);
	Eigen::VectorXd dY = 0.5*(b + c);
#pragma omp parallel for num_threads(24)
	for (int i = 0; i < numF; ++i) {
		//vectors of size 6
		//svd derivatives
		Vector6d dSi = Dsd[0].col(i);
		Vector6d dsi = Dsd[1].col(i);
		//cones constant coefficients (cone = |Ax|, A is a coefficient)
		Vector6d a1i = a1d.col(i);
		Vector6d a2i = a2d.col(i);
		Vector6d b1i = b1d.col(i);
		Vector6d b2i = b2d.col(i);
		Hi[i] = Area(i)*ComputeConvexConcaveFaceHessian(
			a1i, a2i, b1i, b2i,
			aY(i), bY(i), cY(i), dY(i),
			dSi, dsi,
			gradfS(i), gradfs(i),
			HS(i), Hs(i));
		// 		Hi[i].setIdentity();
		int index2 = i * 21;
		for (int a = 0; a < 6; ++a)
		{
			for (int b = 0; b <= a; ++b)
			{
				SS[index2++] = Hi[i](a, b);
			}
		}
	}
	for (int i = 0; i < F.rows(); i++)
	{
		int base = 21 * i;
		SS[base] += 1e-6;
		SS[base + 2] += 1e-6;
		SS[base + 5] += 1e-6;
		SS[base + 9] += 1e-6;
		SS[base + 14] += 1e-6;
		SS[base + 20] += 1e-6;
	}
}

bool ObjectiveAreaPreserving::updateJ(const VectorXd& X)
{
	Eigen::Map<const MatrixX2d> x(X.data(), X.size() / 2, 2);
	// 	a = D1*U;
	// 	b = D2*U;
	// 	c = D1*V;
	// 	d = D2*V;
	for (int i = 0; i < F.rows(); i++)
	{
		// For each face Fi = (p0,p1,p2)
		//		X1i = [ p0.x , p1.x , p2.x ]
		//		X2i = [ p0.y , p1.y , p2.y ]
		Vector3d X1i, X2i;
		X1i << x(F(i, 0), 0), x(F(i, 1), 0), x(F(i, 2), 0);
		X2i << x(F(i, 0), 1), x(F(i, 1), 1), x(F(i, 2), 1);

		// a = Dx.trans * X
		a(i) = D1d.col(i).transpose()*X1i;
		// b = Dy.trans * X
		b(i) = D2d.col(i).transpose()*X1i;
		// c = Dx.trans * Y
		c(i) = D1d.col(i).transpose()*X2i;
		// d = DY.trans * Y
		d(i) = D2d.col(i).transpose()*X2i;

		// detj_1 = det(J) - 1
		double detj_1 = (a(i) * d(i) - b(i) * c(i)) - 1;
		// left = Dx * Y.trans * Dy - Dy * Y.trans * Dx
		Vector3d left = (D1d.col(i) * X2i.transpose() * D2d.col(i)) - (D2d.col(i) * X2i.transpose() * D1d.col(i));
		left *= detj_1;
		// right = Dy * X.trans * Dx - Dx * X.trans * Dy
		Vector3d right = (D2d.col(i) * X1i.transpose() * D1d.col(i)) - (D1d.col(i) * X1i.transpose() * D2d.col(i));
		right *= detj_1;

		grad.row(i) << left(0), left(1), left(2), right(0), right(1), right(2);
	}
	detJ = a.cwiseProduct(d) - b.cwiseProduct(c);
	alpha.col(0) = 0.5*(a + d);   alpha.col(1) = 0.5*(c - b);
	beta.col(0) = 0.5*(a - d);    beta.col(1) = 0.5*(c + b);

	return ((detJ.array() < 0).any());
}

void ObjectiveAreaPreserving::UpdateSSVDFunction()
{

#pragma omp parallel for num_threads(24)
	for (int i = 0; i < a.size(); i++)
	{
		Eigen::Matrix2d A;
		Matrix2d U, S, V;
		A << a[i], b[i], c[i], d[i];
		Utils::SSVD2x2(A, U, S, V);
		u.row(i) << U(0), U(1), U(2), U(3);
		v.row(i) << V(0), V(1), V(2), V(3);
		s.row(i) << S(0), S(3);
	}
}

void ObjectiveAreaPreserving::ComputeDenseSSVDDerivatives()
{
	//different columns belong to diferent faces
	Eigen::MatrixXd B(D1d*v.col(0).asDiagonal() + D2d * v.col(1).asDiagonal());
	Eigen::MatrixXd C(D1d*v.col(2).asDiagonal() + D2d * v.col(3).asDiagonal());

	Eigen::MatrixXd t1 = B * u.col(0).asDiagonal();
	Eigen::MatrixXd t2 = B * u.col(1).asDiagonal();
	Dsd[0].topRows(t1.rows()) = t1;
	Dsd[0].bottomRows(t1.rows()) = t2;
	t1 = C * u.col(2).asDiagonal();
	t2 = C * u.col(3).asDiagonal();
	Dsd[1].topRows(t1.rows()) = t1;
	Dsd[1].bottomRows(t1.rows()) = t2;
}

inline Matrix6d ObjectiveAreaPreserving::ComputeFaceConeHessian(const Vector6d& A1, const Vector6d& A2, double a1x, double a2x)
{
	double f2 = a1x * a1x + a2x * a2x;
	double invf = 1.0 / sqrt(f2);
	double invf3 = invf * invf*invf;

	Matrix6d A1A1t = A1 * A1.transpose();
	Matrix6d A2A2t = A2 * A2.transpose();
	Matrix6d A1A2t = A1 * A2.transpose();
	Matrix6d A2A1t = A1A2t.transpose();


	double a2 = a1x * a1x;
	double b2 = a2x * a2x;
	double ab = a1x * a2x;

	return  (invf - invf3 * a2) * A1A1t + (invf - invf3 * b2) * A2A2t - invf3 * ab*(A1A2t + A2A1t);
}

inline Matrix6d ObjectiveAreaPreserving::ComputeConvexConcaveFaceHessian(const Vector6d& a1, const Vector6d& a2, const Vector6d& b1, const Vector6d& b2, double aY, double bY, double cY, double dY, const Vector6d& dSi, const Vector6d& dsi, double gradfS, double gradfs, double HS, double Hs)
{
	//no multiplying by area in this function

	Matrix6d H = HS * dSi*dSi.transpose() + Hs * dsi*dsi.transpose(); //generalized gauss newton
	double walpha = gradfS + gradfs;
	if (walpha > 0)
		H += walpha * ComputeFaceConeHessian(a1, a2, aY, bY);

	double wbeta = gradfS - gradfs;
	if (wbeta > 1e-7)
		H += wbeta * ComputeFaceConeHessian(b1, b2, cY, dY);
	return H;
}

void ObjectiveAreaPreserving::prepare_hessian()
{
	int n = numV;

	II.clear();
	JJ.clear();
	auto PushPair = [&](int i, int j) { if (i > j) swap(i, j); II.push_back(i); JJ.push_back(j); };
	for (int i = 0; i < F.rows(); ++i)
	{
		// for every face there is a 6x6 local hessian
		// we only need the 21 values contained in the upper
		// triangle. they are access and also put into the
		// big hessian in column order.


		// 		// base indices
		// 		int uhbr = 3 * i; //upper_half_base_row
		// 		int lhbr = 3 * i + n; // lower_half_base_row 
		// 		int lhbc = 3 * i; // left_half_base_col 
		// 		int rhbc = 3 * i + n; // right_half_base_col 
		//uhbr== fi(0)
		Vector3i Fi = F.row(i);
		// first column
		PushPair(Fi(0), Fi(0));

		// second column
		PushPair(Fi(0), Fi(1));
		PushPair(Fi(1), Fi(1));

		// third column
		PushPair(Fi(0), Fi(2));
		PushPair(Fi(1), Fi(2));
		PushPair(Fi(2), Fi(2));

		// fourth column
		PushPair(Fi(0), Fi(0) + n);
		PushPair(Fi(1), Fi(0) + n);
		PushPair(Fi(2), Fi(0) + n);
		PushPair(Fi(0) + n, Fi(0) + n);

		// fifth column
		PushPair(Fi(0), Fi(1) + n);
		PushPair(Fi(1), Fi(1) + n);
		PushPair(Fi(2), Fi(1) + n);
		PushPair(Fi(0) + n, Fi(1) + n);
		PushPair(Fi(1) + n, Fi(1) + n);

		// sixth column
		PushPair(Fi(0), Fi(2) + n);
		PushPair(Fi(1), Fi(2) + n);
		PushPair(Fi(2), Fi(2) + n);
		PushPair(Fi(0) + n, Fi(2) + n);
		PushPair(Fi(1) + n, Fi(2) + n);
		PushPair(Fi(2) + n, Fi(2) + n);
	}
	SS = vector<double>(II.size(), 0.);
}