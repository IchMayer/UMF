#include "FEM.h"

//Константы точности
#define EPS		1e-14	//Погрешность невязки
#define DELTA	1e-14	//Погрешность разности шага решений
#define MAXITER 100000		//Максимальное количество итераций на каждый метод

//X и T
//#define X		{0, 0.25, 0.5, 0.75, 1}
#define T		{0, 0.1, 0.2, 0.3, 0.4}

//Значение констант
#define lambda	1		//Лябда

//Значение функций
#define Sigma_Function	1 //Функция сигма от x t и dudx

//#define	U				pow(x,10) - pow(x,9)	//Искомая функция функция от x и t
//
//#define DuDx			10 * pow(x,9) - 9 * pow(x,8)	//Первая производная функция от x и t
//#define	DivGrad			90 * pow(x,8) - 72 * pow(x,7)	//DivGrad функция от x и t
//
//#define DuDt			0

#pragma region Grid

//#define X		{0, 0.25,0.5, 0.75, 1}
//#define X		{0, 0.125, 0.25, 0.375,0.5,0.675, 0.75,0.875, 1}
//#define X		{0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.675, 0.7375, 0.75, 0.8125, 0.875, 0.9375, 1}

//#define HX2
//#define HX4

//#define X		{0, 0.8, 1}
//#define X		{0, 0.5333333333333333, 0.8, 0.9333333333333333, 1}
//#define X		{0, 0.3555555555555555, 0.5333333333333333, 0.711111111111111, 0.8, 0.888888888888889, 0.9333333333333333, 0.9777777777777777, 1}

//#define T		{0, 0.1, 0.2, 0.3, 0.4}
//#define T		{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4}
//#define T		{0, 0.025, 0.05, 0.75, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4}

//#define H2
//#define H4

#pragma endregion

#pragma region TESTS!
//#define testX
//#define testXX
//#define testXXX
#define testXXXX
//#define testT
//#define testTT
//#define test4XXT
//#define testXT
//#define test3XXXX

#ifdef testX

#define	U				x	//Искомая функция функция от x и t

#define DuDx			1	//Первая производная функция от x и t
#define	DivGrad			0	//DivGrad функция от x и t

#define DuDt			0

#endif

#ifdef testXX

#define	U				x*x	//Искомая функция функция от x и t

#define DuDx			2*x	//Первая производная функция от x и t
#define	DivGrad			2	//DivGrad функция от x и t

#define DuDt			0

#endif

#ifdef testXXX

#define	U				x*x*x	//Искомая функция функция от x и t

#define DuDx			3*x*x	//Первая производная функция от x и t
#define	DivGrad			6*x	//DivGrad функция от x и t

#define DuDt			0

#endif

#ifdef testXXXX

#define	U				x*x*x*x	//Искомая функция функция от x и t

#define DuDx			4*x*x*x	//Первая производная функция от x и t
#define	DivGrad			12*x*x	//DivGrad функция от x и t

#define DuDt			0

#endif

#ifdef testT
#define	U				t	//Искомая функция функция от x и t

#define DuDx			0	//Первая производная функция от x и t
#define	DivGrad			0	//DivGrad функция от x и t

#define DuDt			1 
#endif

#ifdef testTT
#define	U				t*t	//Искомая функция функция от x и t

#define DuDx			0	//Первая производная функция от x и t
#define	DivGrad			0	//DivGrad функция от x и t

#define DuDt			2*t 
#endif

#ifdef testXT
#define	U				x*t	//Искомая функция функция от x и t

#define DuDx			t	//Первая производная функция от x и t
#define	DivGrad			0	//DivGrad функция от x и t

#define DuDt			x 
#endif

#ifdef test4XXT
#define	U				4*x*x*t	//Искомая функция функция от x и t

#define DuDx			8*x*t	//Первая производная функция от x и t
#define	DivGrad			8*t	//DivGrad функция от x и t

#define DuDt			4*x*x 
#endif

#ifdef test3XXXX

#define	U				3*x*x*x*x	//Искомая функция функция от x и t

#define DuDx			12*x*x*x	//Первая производная функция от x и t
#define	DivGrad			36*x*x	//DivGrad функция от x и t

#define DuDt			0 
#endif


#pragma endregion

int main()
{
	FEM m;

	vector<double> x;
	int laksd = 128;
	double yuaysd = 1.0 / laksd;
	for (size_t i = 0; i < laksd; i++)
		x.push_back(i*yuaysd);
	x.push_back(1);

	#pragma region Функции

	//Искомая функция
	function<double(double, double)> u = [](double x, double t) {return U; };
	
	//Првая производная
	function<double(double, double)> dudx = [](double x, double t) {return DuDx; };
	
	//Првая производная по времени
	function<double(double, double)> dudt = [](double x, double t) {return DuDt; };

	//Сумма частных вторых производных
	function<double(double, double)> divGrad = [](double x, double t) {return DivGrad; };

	#pragma endregion
	#pragma region Сигма
	
	function<double(double, double, double)> sigma = [](double x, double t, double dudx) {return Sigma_Function; };
	m.SetSigma(sigma);

	#pragma endregion
	#pragma region Задаем базисы

	FEM::Basis basis;
	basis.type = basis.Lagrange;
	basis.order = 2;

	m.SetBasis(basis);

	#pragma endregion
	#pragma region Задаем константы

	double NX = x.size();

	m.SetLambda(lambda);
	m.eps = EPS;
	m.maxIter = MAXITER;

	#pragma endregion
	#pragma region Задаем краевые условия
	
	function<double(double)> f1 = [&](double t) { return u(x[0], t); };
	function<double(double)> f2 = [&](double t) { return u(x[NX-1], t); };
	FEM::Border b1(f1, FEM::Border::First), b2(f2, FEM::Border::First);

	m.AddBorder(b1);
	m.AddBorder(b2);
	
	#pragma endregion
	#pragma region Задаем форму

	m.SetX(x);
	m.SetT(T);

	#pragma endregion
	#pragma region Задаем начальное значение
	vector<double> q(2 * NX - 1);
	double dx = (x[NX - 1] - x[0]) / (NX - 1);
	for (size_t i = 0; i < NX-1; i++)
	{
		q[2 * i] = u(x[i], vector<double>T[0]);
		q[2 * i + 1] = u((x[i + 1] + x[i]) / 2, vector<double>T[0]);
	}
	q[(NX - 1) * 2] = u(x[NX-1], vector<double>T[0]);

	m.SetQ0(q);
	#pragma endregion
	#pragma region Задаем правую часть уравнения

	function<double(double, double)> F = [&](double x, double t) {return -1 * lambda * divGrad(x, t) + sigma(x, t, dudx(x, t)) * dudt(x, t);};
	m.SetF(F);

	#pragma endregion

	m.Start();

	vector<vector<double>>	Result;
	m.GetResult(Result);

	#pragma region Вывод

#pragma region Определение выводимых строк
	int h = 1, hx = 1;
#ifdef H2
	h = 2;
#endif
#ifdef H4
	h = 4;
#endif
#ifdef HX2
	hx = 2;
#endif
#ifdef HX4
	hx = 4;
#endif

#pragma endregion

	hx = laksd/4;

	cout.precision(11);
	cout << "Значение" << endl;
	cout << setw(6) << "T\\X";

	for (size_t j = 0; j < x.size(); j++)
		if (j % hx == 0)
			cout << setw(18) << x[j];

	cout << "  Итер. в ит. проц."<< endl << endl;

	cout << setw(6) << vector<double>T[0];
	for (size_t j = 0; j < x.size(); j++)
		if (j % hx == 0)
			cout << setw(18) << q[2 *j];

	cout << setw(6) << 0 << endl;

	for (size_t i = 0, i1 = 1; i < Result.size(); i++, i1++)
	{
		if (i1 % h == 0)
		{
			cout << setw(6) << vector<double>T[i + 1];
			for (size_t j = 0; j < x.size(); j++)
				if (j % hx == 0)
					cout << setw(18) << Result[i][2 * j];
			cout << setw(6) << m.It[i] << endl;
		}
	}

	cout << endl << "Погрешность" << endl;

	cout << setw(6) << "T\\X";
	for (size_t j = 0; j < x.size(); j++)
		if (j % hx ==0)
			cout << setw(18) << x[j];

	cout << endl << endl;;

	cout << setw(6) << vector<double>T[0];
	for (size_t j = 0; j < x.size(); j++)
		if (j % hx == 0)
			cout << setw(18) << 0;

	cout << endl;

	for (size_t i = 0, i1 = 1; i < Result.size(); i++, i1++)
	{
		if (i1 % h == 0)
		{
			cout << setw(6) << vector<double>T[i1];
			for (size_t j = 0; j < x.size(); j++)
				if(j % hx == 0)
					cout << setw(18) << u(x[j], vector<double>T[i1]) - Result[i][2 * j];
			cout << endl;
		}
	}
	cout << endl;
	
	//Невязка
	double snev = 0;

	for (size_t i = 0, i1 = 1; i < Result.size(); i++, i1++)
		if (i1 % h == 0)
			for (size_t j = 0; j < x.size(); j++)
				if (j % hx == 0)
					snev += pow(u(x[j], vector<double>T[i1]) - Result[i][2 * j] , 2);

	cout << "Норма вектора невязки: " << sqrt(snev) <<endl;

	#pragma endregion

	system("pause");
	return 0;
}