#include"Includes.h"
#include"Model.h"
//#include <math.h>

#define lamda 1
#define gamma 2
#define betta 2

double lastsum;

void checkstepgrid(int n)
{
	Model a;
	int N = 4;
	int yyy = (n+1) / N;
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

	double hLT = 0.5;
	double hJT = 0.5;

	a.SetKGrid(hLT, hJT);
	a.CreateGrid(n, n);


	a.SetGamma(gamma);
	a.SetLamda(lamda);

	//function < double(double, double) > f = [](double y, double x) { return -1 * lamda * (12 * x * x + 12 * y * y) + gamma * (x * x * x * x + y * y*y * y); };	//-lamda div grad(U) + gamma * U = f
	//function < double(double, double) > u = [](double y, double x) { return x * x * x * x + y * y * y * y; };							// U
	//function < double(double, double) > q = [](double y, double x) { return -1 * lamda * (2 * x); };			//  -lamda dU/dn
	//function < double(double, double) > q2 = [](double y, double x) { return -1 * lamda * (2 * y); };
	//function < double(double, double) > ub = [](double y, double x) { return lamda * 4 * y * y * y + betta * (pow(x, 4) + pow(y, 4)); };
	//function < double(double, double) > ub2 = [](double y, double x) { return lamda * 4 * x * x * x + betta * (pow(x, 4) + pow(y, 4)); };

	//function < double(double, double) > f = [](double y, double x) { return -1 * lamda * (4) + gamma * (x * x + y * y); };	//-lamda div grad(U) + gamma * U = f
	//function < double(double, double) > u = [](double y, double x) { return x * x + y * y; };							// U
	//function < double(double, double) > q = [](double y, double x) { return -1 * lamda * (2 * x); };			//  -lamda dU/dn
	//function < double(double, double) > q2 = [](double y, double x) { return -1 * lamda * (2 * y); };
	//function < double(double, double) > ub = [](double y, double x) { return lamda * 2 * y + betta * (pow(x, 2) + pow(y, 2)); };
	//function < double(double, double) > ub2 = [](double y, double x) { return lamda * 2 * x + betta * (pow(x, 2) + pow(y, 2)); };

	function < double(double, double) > f = [](double y, double x) { return -1 * lamda * (2) + gamma * (x * x); };	//-lamda div grad(U) + gamma * U = f
	function < double(double, double) > u = [](double y, double x) { return x * x; };							// U
	function < double(double, double) > q = [](double y, double x) { return -1 * lamda * (2 * x); };			//  -lamda dU/dn
	function < double(double, double) > q2 = [](double y, double x) { return 0; };
	function < double(double, double) > ub = [](double y, double x) { return lamda * 0 + betta * (pow(x, 2) + pow(y, 2)); };
	function < double(double, double) > ub2 = [](double y, double x) { return lamda * 2 * x + betta * (pow(x, 2) + pow(y, 2)); };

	//function < double(double, double) > f = [](double y, double x) { return -1 * lamda * (exp(x) + exp(y)) + gamma * (exp(x) + exp(y)); };	//-lamda div grad(U) + gamma * U = f
	//function < double(double, double) > u = [](double y, double x) { return exp(x) + exp(y); };							// U
	//function < double(double, double) > q = [](double y, double x) { return -1 * lamda * (exp(x)); };			//  -lamda dU/dn
	//function < double(double, double) > q2 = [](double y, double x) { return -1 * lamda * (exp(y)); };
	//function < double(double, double) > ub = [](double y, double x) { return lamda * exp(y) + betta * (exp(x) + exp(y)); };
	//function < double(double, double) > ub2 = [](double y, double x) { return lamda * exp(x) + betta * (exp(x) + exp(y)); };

	Model::Border b1;

	b1.type = b1.First;
	b1.borderF = u;

	Model::Border b2;
	Model::Border b22;
	Model::Border b3;
	Model::Border b4;

	b3.type = b3.Third;
	b4.type = b4.Third;
	b3.borderF = ub2;
	b4.borderF = ub;
	b3.b = betta;
	b4.b = betta;

	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);

	a.SetF(f);

	a.FDM();

	vector<vector<double>> U = a.GetResult();

	double hHI = a.GethHI();
	double hHJ = a.GethHJ();

	vector<vector<double>> UR(U.size());

	for (size_t i = 0; i < U.size(); i++)
	{
		UR[i].resize(U.size());
		for (size_t j = 0; j < U[i].size(); j++)
		{
			double sum1 = shape.minX;
			double sum2 = shape.minY;

			double kI = hHI;
			double kJ = hHJ;

			for (size_t k1 = 0; k1 < i; k1++, kI *= hLT)
				sum1 += kI;
			for (size_t k1 = 0; k1 < j; k1++, kJ *= hJT)
				sum2 += kJ;

			UR[i][j] = u(sum1, sum2);
		}
	}

	double max = 0;
	double sum = 0;
	double d;
	double sumN1 = 0;
	for (size_t i = 0, i1 = n; i < n + 1; i++, i1--)
	{
		for (size_t j = 0; j < n + 1; j++)
		{
			if (i % yyy == 0 && j % yyy == 0)
			{
				d = abs(U[i][j] ? U[i][j] - UR[i1][j] : 0);
				sumN1 += pow(d, 2);
				sum += d;
				if (d > max)
					max = d;
			}		
		}
	}

	cout << n + 1<< " \t" << /*max <<*/ " \t"<< sqrt(sumN1) << endl;
	
	//if (lastsum)
	//	cout << n + 1<< "\t\t " << log2(lastsum / sqrt(sumN1)) << endl;
	//lastsum = sqrt(sumN1);
}

int main()
{
	lastsum = 0;
	for (size_t i = 2; i < 8; i++)
	{
		checkstepgrid(pow(2, i) - 1);
	}

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
	a.SetKGrid(0.5, 0.5);
	a.CreateGrid(15, 15);
	a.SetGamma(gamma);
	a.SetLamda(lamda);
	
	//function < double(double, double) > f = [](double y, double x) { return -4; };
	//function < double(double, double) > u = [](double y, double x) { return x * x + y * y; };
	//function < double(double, double) > q = [](double y, double x) { return -2 * x; };
	//function < double(double, double) > q2 = [](double y, double x) { return -2 * y; };
	//function < double(double, double) > ub = [](double y, double x) { return 2 * y + 2 *(x * x + y * y); };
	//function < double(double, double) > ub2 = [](double y, double x) { return 2 * x + 2 *(x * x + y * y); };
	
	function < double(double, double) > f = [](double y, double x) { return -1 * lamda * (4) + gamma * (x * x + y * y); };	//-lamda div grad(U) + gamma * U = f
	function < double(double, double) > u = [](double y, double x) { return x * x + y * y; };							// U
	function < double(double, double) > q = [](double y, double x) { return -1 * lamda * (2 * x); };			//  -lamda dU/dn
	function < double(double, double) > q2 = [](double y, double x) { return -1 * lamda * (2 * y); };
	function < double(double, double) > ub = [](double y, double x) { return lamda * 2 * y + betta * (pow(x, 2) + pow(y, 2)); };
	function < double(double, double) > ub2 = [](double y, double x) { return lamda * 2 * x + betta * (pow(x, 2) + pow(y, 2)); };


	Model::Border b1;

	b1.type = b1.First;
	b1.borderF = u;

	Model::Border b2;
	Model::Border b22;
	Model::Border b3;
	Model::Border b4;

	b2.type = b2.Second;
	b22.type = b22.Second;
	b2.borderF = q;
	b22.borderF = q2;

	b3.type = b3.Third;
	b4.type = b4.Third;
	b3.borderF = ub2;
	b4.borderF = ub;
	b3.b = betta;
	b4.b = betta;

	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);
	a.AddBorder(b1);

	a.SetF(f);

	a.FDM();

	vector<vector<double>> U = a.GetResult();

	for (size_t i = 0; i < U.size(); i++)
	{
		for (size_t j = 0; j < U[i].size(); j++)
		{
			cout.width(15);
			cout << U[i][j];
		}
		cout << endl;
	}

	cout << endl;
	cout << endl;
	cout << endl;

	double hLT = 0.5;
	double hJT = 0.5;

	double hHI = a.GethHI();
	double hHJ = a.GethHJ();

	vector<vector<double>> UR(U.size());

	for (size_t i = 0; i < U.size(); i++)
	{
		UR[i].resize(U.size());
		for (size_t j = 0; j < U[i].size(); j++)
		{
			double sum1 = shape.minX;
			double sum2 = shape.minY;

			double kI = hHI;
			double kJ = hHJ;

			for (size_t k1 = 0; k1 < i; k1++, kI *= hLT)
				sum1 += kI;
			for (size_t k1 = 0; k1 < j; k1++, kJ *= hJT)
				sum2 += kJ;

			UR[i][j] = u(sum1, sum2);
		}
	}

	for (size_t i = 0, i1 = U.size() - 1; i < U.size(); i++, i1--)
	{
		for (size_t j = 0; j < U.size(); j++)
		{
			cout.width(15);
			cout << (U[i][j] ? U[i][j] - UR[i1][j] : 0) << " ";
		}
		cout << endl;
	}

	int i;
	cin >> i;

	return 0;
}