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
	shape.type = shape.Rectangle;

	a.SetShape(shape);
	a.CreateGrid(4, 4);
	
	function < double(double, double) > f = [](double x, double y) { return 2 * x + 2 * y; };
	function < double(double, double) > u = [](double x, double y) { return x + y; };

	Model::Border b1;

	b1.type = b1.First;
	b1.borderF = u;

	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);

	a.SetF(f);

	a.FDM();

	return 0;
}