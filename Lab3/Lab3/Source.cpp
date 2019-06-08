#include "FEM.h"

//Константы точности
#define EPS		1e-14	//Погрешность невязки
#define DELTA	1e-14	//Погрешность разности шага решений
#define MAXITER 100000		//Максимальное количество итераций на каждый метод

//X и T
#define X		{0, 0.25, 0.5, 0.75, 1}
#define Y		{0, 0.25, 0.5, 0.75, 1}
#define T		{0, 0.1, 0.2, 0.3, 0.4}

//Значение констант
#define lambda	1		//Лябда
#define sigma	1 //Функция сигма от x t и dudx
#define omega 1
#define hi 1

//U функция от x y, t
#define Us x
#define Uc x

#define DivGrads 0
#define DivGradc 0

int main()
{
	FEM m;

	vector<double> x = X;
	vector<double> y = Y;

	#pragma region Функции

	//Искомая функция
	function<double(double, double)> us = [](double x, double y) { return Us; };
	function<double(double, double)> uc = [](double x, double y) { return Uc; };

	function<double(double, double)> divgradS = [](double x, double y) {return DivGrads; };
	function<double(double, double)> divgradC = [](double x, double y) {return DivGradc; };

	#pragma endregion
	#pragma region Задаем базисы

	FEM::Basis basis;
	basis.type = basis.Lagrange;
	basis.order = 1;

	m.SetBasis(basis);

	#pragma endregion
	#pragma region Задаем константы

	int NX = x.size();
	int NY = y.size();

	m.SetLambda(lambda);
	m.SetSigma(sigma);
	m.SetOmega(omega);
	m.eps = EPS;
	m.maxIter = MAXITER;

	#pragma endregion
	#pragma region Задаем краевые условия
	
	function<double(double, double)> f1s = [&](double y, double t) { return us(x[0], y); };
	function<double(double, double)> f1c = [&](double y, double t) { return uc(x[0], y); };
	function<double(double, double)> f2s = [&](double x, double t) { return us(x, y[0]); };
	function<double(double, double)> f2c = [&](double x, double t) { return uc(x, y[0]); };
	function<double(double, double)> f3s = [&](double y, double t) { return us(x[NX-1], y); };
	function<double(double, double)> f3c = [&](double y, double t) { return uc(x[NX-1], y); };
	function<double(double, double)> f4s = [&](double x, double t) { return us(x, y[NY-1]); };
	function<double(double, double)> f4c = [&](double x, double t) { return uc(x, y[NY-1]); };
	
	FEM::Border b1(f1s, f1c, FEM::Border::First),
				b2(f2s, f2c, FEM::Border::First),
				b3(f3s, f3c, FEM::Border::First),
				b4(f4s, f4c, FEM::Border::First);

	m.AddBorder(b1);
	m.AddBorder(b2);
	m.AddBorder(b3);
	m.AddBorder(b4);
	
	#pragma endregion
	#pragma region Задаем форму

	m.SetX(x);
	m.SetY(y);
	m.SetT(T);

	#pragma endregion
	#pragma region Задаем начальное значение ЕЩЕ НЕТ

	#pragma endregion
	#pragma region Задаем правую часть уравнения

	
	function<double(double, double, double)> Fs = [&](double x, double y) {return -1 * divgradS(x, y) - omega * sigma * uc(x, y) - pow(omega,2) * hi * us(x,y)};
	function<double(double, double, double)> Fc = [&](double x, double y) {return -1 * divgradC(x, y) + omega * sigma * us(x, y) - pow(omega,2) * hi * uc(x, y);};
	
	m.SetF(Fc, Fs);

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

	//hx = laksd/4;

	//cout.precision(11);
	//cout << "Значение" << endl;
	//cout << setw(6) << "T\\X";

	//for (size_t j = 0; j < x.size(); j++)
	//	if (j % hx == 0)
	//		cout << setw(18) << x[j];

	//cout << "  Итер. в ит. проц."<< endl << endl;

	//cout << setw(6) << vector<double>T[0];
	//for (size_t j = 0; j < x.size(); j++)
	//	if (j % hx == 0)
	//		cout << setw(18) << q[2 *j];

	//cout << setw(6) << 0 << endl;

	//for (size_t i = 0, i1 = 1; i < Result.size(); i++, i1++)
	//{
	//	if (i1 % h == 0)
	//	{
	//		cout << setw(6) << vector<double>T[i + 1];
	//		for (size_t j = 0; j < x.size(); j++)
	//			if (j % hx == 0)
	//				cout << setw(18) << Result[i][2 * j];
	//		cout << setw(6) << m.It[i] << endl;
	//	}
	//}

	//cout << endl << "Погрешность" << endl;

	//cout << setw(6) << "T\\X";
	//for (size_t j = 0; j < x.size(); j++)
	//	if (j % hx ==0)
	//		cout << setw(18) << x[j];

	//cout << endl << endl;;

	//cout << setw(6) << vector<double>T[0];
	//for (size_t j = 0; j < x.size(); j++)
	//	if (j % hx == 0)
	//		cout << setw(18) << 0;

	//cout << endl;

	//for (size_t i = 0, i1 = 1; i < Result.size(); i++, i1++)
	//{
	//	if (i1 % h == 0)
	//	{
	//		cout << setw(6) << vector<double>T[i1];
	//		for (size_t j = 0; j < x.size(); j++)
	//			if(j % hx == 0)
	//				cout << setw(18) << u(x[j], vector<double>T[i1]) - Result[i][2 * j];
	//		cout << endl;
	//	}
	//}
	//cout << endl;
	//
	////Невязка
	//double snev = 0;

	//for (size_t i = 0, i1 = 1; i < Result.size(); i++, i1++)
	//	if (i1 % h == 0)
	//		for (size_t j = 0; j < x.size(); j++)
	//			if (j % hx == 0)
	//				snev += pow(u(x[j], vector<double>T[i1]) - Result[i][2 * j] , 2);

	//cout << "Норма вектора невязки: " << sqrt(snev) <<endl;

	#pragma endregion

	system("pause");
	return 0;
}