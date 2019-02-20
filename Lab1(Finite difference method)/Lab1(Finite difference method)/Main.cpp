#include"Includes.h"
#include"Model.h"

int main()
{
	Model a;
	Model::Shape shape;
	shape.maxX = 1;
	shape.maxY = 1;
	shape.minX = 0;
	shape.minY = 0;
	shape.aXT = 0.75;
	shape.iXT = 0.25;
	shape.YT = 0.5;
	shape.type = shape.T;

	a.SetShape(shape);
	a.CreateGrid(4, 4);
	
	function < double(double, double) > f = [](double y, double x) { return 0; };
	function < double(double, double) > u = [](double y, double x) { return 2 * x + 2 * y; };
	function < double(double, double) > q = [](double y, double x) { return -2; };
	function < double(double, double) > ub = [](double y, double x) { return 2 + 2 *(2 * x + 2 * y); };
	
	Model::Border b1;

	b1.type = b1.First;
	b1.borderF = u;

	Model::Border b2;
	Model::Border b3;

	b2.type = b2.Second;
	b2.borderF = q;

	b3.type = b3.Third;
	b3.borderF = ub;
	b3.b = 2;

	a.AddBorder(b3);
	a.AddBorder(b3);
	a.AddBorder(b3);
	a.AddBorder(b3);
	a.AddBorder(b3);
	a.AddBorder(b3);
	a.AddBorder(b3);
	a.AddBorder(b3);

	a.SetF(f);

	a.FDM();

	vector<vector<double>> U = a.GetResult();

	for (size_t i = 0; i < U.size(); i++)
	{
		for (size_t j = 0; j < U[i].size(); j++)
		{
			cout.width(10);
			cout << U[i][j];
		}
		cout << endl;
	}

	int i;
	cin >> i;

	return 0;
}