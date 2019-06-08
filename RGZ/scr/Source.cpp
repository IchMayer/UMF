#include "FEM.h"

//Константы точности
#define EPS		1e-14	//Погрешность невязки
#define DELTA	1e-14	//Погрешность разности шага решений
#define MAXITER 100000		//Максимальное количество итераций на каждый метод

//X и T
#define X		{1, 3, 5, 7}
#define Y		{1, 2, 3}
#define T		{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}

//Значение констант
#define lambda	3		//Лябда
#define sigma	1 //Функция сигма от x t и dudx
#define hi 0.1

#define file "input.txt"

//3 краевое
//#define betta 1

//#define U x * x * t * t * t + y * y * y * t * t
//
//#define DUDT 3 * x * x * t * t + 2 * y * y * y * t
//#define D2UDT2 6 * x * x * t + 2 * y * y * y
//
//#define DUDX 2 * x * t * t * t
//#define DUDY 3 * y * y * t * t
//
//#define DivGrad 2 * t * t * t + 6 * y * t * t 

//#define Test1
//#define Test2
//#define Test3
//#define Test4
#define Test5
#ifdef Test1

#define U 3 * pow(x, 2) + 2 * pow(y, 3)
#define DUDT 0
#define D2UDT2 0

#define DUDX 6 * x 
#define DUDY 6 * pow(y, 2)

#define DivGrad 6 + 12 * y

#define Bord {1, 1, 1, 1}

#define betta 1
#endif

#ifdef Test2

#define U pow(y, 4)
#define DUDT 0
#define D2UDT2 0

#define DUDX 0 
#define DUDY 4 * pow(y, 3)

#define DivGrad 12 * pow(y, 2)

#define Bord {2, 1, 2, 2}

#define betta 1
#endif

#ifdef Test3

#define U 3 * pow(x, 2) * pow(t, 3) + 2 * pow(y, 3) * pow(t, 2)
#define DUDT 9 * pow(x, 2) * pow(t, 2) + 4 * pow(y, 3) * t
#define D2UDT2 18 * pow(x, 2) * t + 4 * pow(y, 3)

#define DUDX 6 * x * pow(t, 3) 
#define DUDY 6 * pow(y, 2) * pow(t, 2)

#define DivGrad 6 * pow(t, 3) + 12 * y * pow(t, 2)

#define Bord {1, 1, 1, 1}

#define betta 1
#endif

#ifdef Test4

#define U pow(t, 4)
#define DUDT 4 * pow(t, 3)
#define D2UDT2 12 * pow(t, 2)

#define DUDX 0 
#define DUDY 0

#define DivGrad 0

#define Bord {1, 1, 1, 1}

#define betta 1

#endif

#ifdef Test5

#define U exp(4*x*t)
#define DUDT 4*x*exp(4*x*t)
#define D2UDT2 16*x*x*exp(4*x*t)

#define DUDX 4*t*exp(4*x*t)
#define DUDY 0

#define DivGrad 16*t*t*exp(4*x*t)

#define Bord {1, 3, 3, 3}

#define betta 1
#endif

//#define Print

vector<double> createwline(double w, double a, double b, double n)
{
	vector<double> c;
	double r = (b - a)* (1 - w) / (1 - pow(w, n));
	c.push_back(a);
	for (size_t i = 0; i < n; i++)
	{
		r *= w;
		a += r;
		c.push_back(a);
	}
	return c;
}

vector<double> createline(double a, double b, double n)
{
	vector<double> c;
	double r = (b - a) / n;
	for (size_t i = 0; i < n + 1; i++)
		c.push_back(a + i * r);
	return c;
}

int main()
{
	FEM m;

	int n = 128;
	int n0 = n;
	int nt = 8;

#ifdef  file
	ifstream in(file);

	int sizex;
	int sizey;
	int sizet;

	in >> sizex >> sizey >> sizet;

	double clamba;
	double csigma;
	double chi;

	in >> clamba >> csigma >> chi;

	vector<double> x(sizex);
	vector<double> y(sizey);
	vector<double> t(sizet);

	for (size_t i = 0; i < sizex; i++)
		in >> x[i];

	for (size_t i = 0; i < sizey; i++)
		in >> y[i];

	for (size_t i = 0; i < sizet; i++)
		in >> t[i]

#else
	vector<double> x = createline(0, 1, 10);
	vector<double> y = createline(0, 1, 10);
	vector<double> t = createline(0, 1, n);



#endif //  file

	vector<int> bord = Bord;
	//for (size_t i = 0; i < 100; i++)
	//{
	//	x.push_back(i);
	//	y.push_back(i);
	//}

	#pragma region Функции

	//Искомая функция
	function<double(double, double, double)> u = [](double x, double y, double t) { return U; };
	function<double(double, double, double)> dudt = [](double x, double y, double t) {return DUDT; };
	function<double(double, double, double)> d2udt2 = [](double x, double y, double t) {return D2UDT2; };
	function<double(double, double, double)> divgrad = [](double x, double y, double t) {return DivGrad; };

	function<double(double, double, double)> dudx = [](double x, double y, double t) {return DUDX; };
	function<double(double, double, double)> dudy = [](double x, double y, double t) {return DUDY; };

	#pragma endregion
	#pragma region Задаем базисы

	FEM::Basis basis;
	basis.type = basis.Lagrange;
	basis.order = 2;

	m.SetBasis(basis);

	#pragma endregion
	#pragma region Задаем константы

	int NX = x.size();
	int NY = y.size();

	m.SetLambda(lambda);
	m.SetSigma(sigma);
	m.SetHi(hi);
	m.eps = EPS;
	m.maxIter = MAXITER;

	#pragma endregion
	#pragma region Задаем краевые условия
	
	function<double(double, double, double t)> f1 = [&](double z, double y, double t) { return u(x[0], y, t)	; };
	function<double(double, double, double t)> f2 = [&](double x, double z, double t) { return u(x, y[0], t)	; };
	function<double(double, double, double t)> f3 = [&](double z, double y, double t) { return u(x[NX-1], y, t)	; };
	function<double(double, double, double t)> f4 = [&](double x, double z, double t) { return u(x, y[NY-1], t)	; };
	
	function<double(double, double, double t)> f12 = [&](double z, double y, double t) { return lambda * (-dudx(x[0], y, t)); };
	function<double(double, double, double t)> f22 = [&](double x, double z, double t) { return lambda * (-dudy(x, y[0], t)); };
	function<double(double, double, double t)> f32 = [&](double z, double y, double t) { return lambda * (dudx(x[NX - 1], y, t)); };
	function<double(double, double, double t)> f42 = [&](double x, double z, double t) { return lambda * (dudy(x, y[NY - 1], t)); };

	function<double(double, double, double t)> f13 = [&](double z, double y, double t) { return lambda / betta * (-dudx(x[0], y, t))		+ u(x[0], y, t); };
	function<double(double, double, double t)> f23 = [&](double x, double z, double t) { return lambda / betta * (-dudy(x, y[0], t))		+ u(x, y[0], t); };
	function<double(double, double, double t)> f33 = [&](double z, double y, double t) { return lambda / betta * (dudx(x[NX - 1], y, t))	+ u(x[NX - 1], y, t); };
	function<double(double, double, double t)> f43 = [&](double x, double z, double t) { return lambda / betta * (dudy(x, y[NY - 1], t))	+ u(x, y[NY - 1], t); };

	FEM::Border b1(f1, FEM::Border::First),
				b2(f2, FEM::Border::First),
				b3(f3, FEM::Border::First),
				b4(f4, FEM::Border::First),
				b12(f12, FEM::Border::Second),
				b22(f22, FEM::Border::Second),
				b32(f32, FEM::Border::Second),
				b42(f42, FEM::Border::Second),
				b13(f13, FEM::Border::Third, betta),
				b23(f23, FEM::Border::Third, betta),
				b33(f33, FEM::Border::Third, betta),
				b43(f43, FEM::Border::Third, betta);;

	switch (bord[0])
	{
	case 1: m.AddBorder(b1); break;
	case 2: m.AddBorder(b12); break;
	case 3: m.AddBorder(b13); break;
	default: m.AddBorder(b1); break;
	}

	switch (bord[1])
	{
	case 1: m.AddBorder(b2); break;
	case 2: m.AddBorder(b22); break;
	case 3: m.AddBorder(b23); break;
	default:m.AddBorder(b2); break;
	}

	switch (bord[2])
	{
	case 1: m.AddBorder(b3); break;
	case 2: m.AddBorder(b32); break;
	case 3: m.AddBorder(b33); break;
	default: m.AddBorder(b3); break;
	}

	switch (bord[3])
	{
	case 1: m.AddBorder(b4); break;
	case 2: m.AddBorder(b42); break;
	case 3: m.AddBorder(b43); break;
	default: m.AddBorder(b4); break;
	}
	
	#pragma endregion
	#pragma region Задаем форму

	m.SetX(x);
	m.SetY(y);

	m.SetT(t);
	#pragma endregion
	#pragma region Задаем правую часть уравнения

	function<double(double, double, double)> F = [&](double x, double y, double t) {return lambda * -1 * divgrad(x, y, t) + sigma * dudt(x, y, t) + hi * d2udt2(x, y, t); };

	m.SetF(F);

	#pragma endregion
	#pragma region Задаем начальное значени

	for (size_t k = 0; k < 3; k++)
	{
		vector<double> q;
		for (size_t i = 0, i1 = 0; i < y.size(); i++)
		{
			for (size_t j = 0; j < x.size(); j++, i1++)
			{
				q.push_back(u(x[j], y[i], t[k]));
				if (j + 1 < x.size())
					q.push_back(u((x[j] + x[j + 1]) / 2, y[i], t[k]));
			}

			if(i + 1 < y.size())
			for (size_t j = 0; j < x.size(); j++, i1++)
			{
				q.push_back(u(x[j], (y[i] + y[i + 1]) / 2, t[k]));
				if(j + 1 < x.size())
					q.push_back(u((x[j] + x[j + 1]) / 2, (y[i] + y[i + 1]) / 2, t[k]));
			}
		}
		m.AddSolution(q);
	}

	#pragma endregion



	#pragma region Вывод
	
	m.Start();

	vector<vector<double>> res;

	m.GetResult(res);

	int w = 15;

#ifdef Print
	cout << "Решение" << endl;

	for (size_t k = 0; k < res.size(); k++)
	{
		cout << "t = " << t[k] << endl;
		cout << setw(w) << "y\\x";

		for (size_t j = 0; j < 2 * x.size() - 1; j++)
		{
			if (j % 2)
				cout << setw(w) << (x[j / 2] + x[(j + 1) / 2]) / 2 << " ";
			else
				cout << setw(w) << x[j / 2] << " ";
		}
		cout << endl;

		for (size_t i = 0, i1 = 0; i < 2 * y.size() - 1; i++)
		{
			if (i % 2)
				cout << setw(w) << (y[i / 2] + y[(i + 1) / 2]) / 2 << " ";
			else
				cout << setw(w) << y[i / 2] << " ";

			for (size_t j = 0; j < 2 * x.size() - 1; j++, i1++)
				cout << setw(w) << res[k][i1] << " ";
			cout << endl;
		}
		cout << endl;
	}

	cout << endl << "Погрешность" << endl;

	for (size_t i = 0; i < res.size(); i++)
	{
		cout << "t = " << t[i] << endl;
		cout << setw(w) << "y\\x";

		for (size_t j = 0; j < 2 * x.size() - 1; j++)
		{
			if (j % 2)
				cout << setw(w) << (x[j / 2] + x[(j + 1) / 2]) / 2 << " ";
			else
				cout << setw(w) << x[j / 2] << " ";
		}
		cout << endl;
		for (size_t j = 0, j1 = 0; j < 2 * y.size() - 1; j++)
		{
			{
				if (j % 2)
				{
					cout << setw(w) << (y[j / 2] + y[(j + 1) / 2]) / 2 << " ";
					for (size_t k = 0; k < 2 * x.size() - 1; k++, j1++)
						if (k % 2)
							cout << setw(w) << u((x[k / 2] + x[(k + 1) / 2]) / 2, (y[j / 2] + y[(j + 1) / 2]) / 2, t[i]) - res[i][j1] << " ";
						else
							cout << setw(w) << u(x[k / 2], (y[j / 2] + y[(j + 1) / 2]) / 2, t[i]) - res[i][j1] << " ";
				}
				else
				{
					cout << setw(w) << y[j / 2] << " ";
					for (size_t k = 0; k < 2 * x.size() - 1; k++, j1++)
						if (k % 2)
							cout << setw(w) << u((x[k / 2] + x[(k + 1) / 2]) / 2, y[j / 2], t[i]) - res[i][j1] << " ";
						else
							cout << setw(w) << u(x[k / 2], y[j / 2], t[i]) - res[i][j1] << " ";
				}
				cout << endl;
			}
		}
		normvect.push_back(sqrt(norm));
		norm = 0;
		cout << endl;
	}

#endif // Print

	double norm = 0;
	double normt = 0;
	vector<double> normvect(0);


	for (size_t i = 0; i < res.size(); i++)
	{
		for (size_t j = 0, j1 = 0; j < 2 * y.size() - 1; j++)
		{
			if (j % 2)
			{
				for (size_t k = 0; k < 2 * x.size() - 1; k++, j1++)
					if (k % 2)
					{
						if (j % (n / n0) == 0 && k % (n / n0) == 0)	norm += pow(u((x[k / 2] + x[(k + 1) / 2]) / 2, (y[j / 2] + y[(j + 1) / 2]) / 2, t[i]) - res[i][j1], 2);
					}
					else
					{
						if (j % (n / n0) == 0 && k % (n / n0) == 0) norm += pow(u(x[k / 2], (y[j / 2] + y[(j + 1) / 2]) / 2, t[i]) - res[i][j1], 2);
					}
			}
			else
			{
				for (size_t k = 0; k < 2 * x.size() - 1; k++, j1++)
					if (k % 2)
					{
						if (j % (n / n0) == 0 && k % (n / n0) == 0) norm += pow(u((x[k / 2] + x[(k + 1) / 2]) / 2, y[j / 2], t[i]) - res[i][j1], 2);
					}
					else
					{
						if (j % (n / n0) == 0 && k % (n / n0) == 0) norm += pow(u(x[k / 2], y[j / 2], t[i]) - res[i][j1], 2);
					}
			}
			//cout << endl;
		}
		if (i % (n / nt) == 0) normt += norm;

		normvect.push_back(sqrt(norm));
		norm = 0;
	}

	std::cout << "!!!!!!!!! " << sqrt(normt) << endl;

	std::cout << "Норма вектора погрешности:  " << endl;

	for (size_t i = 0; i < t.size(); i++)
	{
	//	std::cout << "t = " << t[i] << "    " << normvect[i] << endl;
	}

	#pragma endregion

	system("pause");
	return 0;
}