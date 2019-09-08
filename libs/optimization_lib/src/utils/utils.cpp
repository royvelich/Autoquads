#include "utils/utils.h"

void Utils::SSVD2x2(const Eigen::Matrix2d& A, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (A(0) + A(3)) * 0.5;
	double f = (A(0) - A(3)) * 0.5;
	double g = (A(1) + A(2)) * 0.5;
	double h = (A(1) - A(2)) * 0.5;
	double q = sqrt((e * e) + (h * h));
	double r = sqrt((f * f) + (g * g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1) * 0.5;
	double phi = (a2 + a1) * 0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}