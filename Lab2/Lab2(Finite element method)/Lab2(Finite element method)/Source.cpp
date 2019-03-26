#include "FEM.h"

int main()
{
	FEM m;

	FEM::Basis basis;
	basis.type = basis.Lagrange;
	basis.order = 2;

	m.SetBasis(basis);

	double gamma = 0.5;
	double lambda = 1;
	double t0 = 0;
	double nt = 10;
	double ht = 0.1;


	m.SetLambda(lambda);
	m.SetGamma(gamma);
	m.eps = 1e-14;
	m.maxIter = 10;
	
	m.SetHt(0.5);
	
	function<double(double)> f1 = [](double t) { return 0; };
	function<double(double)> f2 = [](double t) { return 1; };
	FEM::Border b1(f1, FEM::Border::First), b2(f2, FEM::Border::First);  

	m.AddBorder(b1);
	m.AddBorder(b2);

	FEM::Shape shape;
	shape.n = 5;
	shape.xMax = 0;
	shape.xMax = 1;
	shape.type = shape.Recrangle;

	m.SetShape(shape);

	function<double(double, double)> F = [&](double x, double t) {return gamma * x; };
	m.SetF(F);

	m.SetT0(t0);
	m.SetHt(ht);
	m.SetNt(nt);

	m.Start();

	return 0;
}