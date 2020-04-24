#include "objective_functions/AllVertexPositions.h"

AllVertexPositions::AllVertexPositions()
{
    name = "All Vertex Positions";
	w = 10000;
}

void AllVertexPositions::init()
{
	if (restShapeV.size() == 0 || restShapeF.size() == 0)
		throw name + " must define members V,F before init()!";
	init_hessian();
}

void AllVertexPositions::updateX(const Eigen::VectorXd& X)
{
	assert(X.rows() == (3 * restShapeV.rows()));
	CurrV = Eigen::Map<const Eigen::MatrixX3d>(X.data(), X.rows() / 3, 3);
}

double AllVertexPositions::value(const bool update)
{
	double E = (CurrV - restShapeV).squaredNorm();
	if (update)
		energy_value = E;
	return E;
}

void AllVertexPositions::gradient(Eigen::VectorXd& g, const bool update)
{
	int n = restShapeV.rows();
	g.conservativeResize(n*3);
	g.setZero();

	Eigen::MatrixX3d diff = CurrV - restShapeV;
	for (int i = 0; i < n; i++) {
		g(i + (0 * n)) = 2 * diff(i, 0); //X-coordinate
		g(i + (1 * n)) = 2 * diff(i, 1); //Y-coordinate
		g(i + (2 * n)) = 2 * diff(i, 2); //Z-coordinate
	}
	if(update)
		gradient_norm = g.norm();
}

void AllVertexPositions::hessian()
{
	int n = restShapeV.rows();
	fill(SS.begin(), SS.end(), 0);
	for (int i = 0; i < n; i++)
	{
		SS[i + (0 * n)] = 2;
		SS[i + (1 * n)] = 2;
		SS[i + (2 * n)] = 2;
	}
}

void AllVertexPositions::init_hessian()
{
	int n = restShapeV.rows();
	II.resize(3 * n);
	JJ.resize(3 * n);
	for (int i = 0; i < 3*n; i++)
	{
		II[i] = i;
		JJ[i] = i;
	}
	SS = std::vector<double>(II.size(), 0.);
}