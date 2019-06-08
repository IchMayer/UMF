#include "FEM.h"

//Константы точности
#define EPS		1e-14	//Погрешность невязки
#define DELTA	1e-14	//Погрешность разности шага решений
#define MAXITER 100000		//Максимальное количество итераций на каждый метод

//X и T
#define X		{0, 0.78539816339, 1.57079632679, 2.35619449019, 3.14159265359}
#define Y		{0, 0.78539816339, 1.57079632679, 2.35619449019, 3.14159265359}

//Значение констант
#define lambda	1		//Лябда
#define sigma	1 //Функция сигма от x t и dudx
#define omega 1
#define hi 1

//U функция от x y, t
//#define Us pow(x,3) * y
//#define Uc pow(y,3) * x
//
//#define DivGrads 6 * x * y
//#define DivGradc 6 * y * x

//#define Us pow(x,3)*y
//#define Uc pow(y,3)*x
//
//#define DivGrads 6 * x * y
//#define DivGradc 6 * y * x


#define Us 0
#define Uc 1

#define DivGrads 0
#define DivGradc 0

//#define Check
//#define CheckLambda
//#define CheckSigma
//#define CheckOmega
//#define CheckHI
//#define CheckEPS
#define Print

int main()
{
	FEM m;

	vector<double> x = X;
	vector<double> y = Y;

	//for (size_t i = 0; i < 100; i++)
	//{
	//	x.push_back(i);
	//	y.push_back(i);
	//}

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
	m.SetHi(hi);
	m.eps = EPS;
	m.maxIter = MAXITER;

	#pragma endregion
	#pragma region Задаем краевые условия
	
	function<double(double, double)> f1s = [&](double z, double y) { return us(x[0], y); };
	function<double(double, double)> f1c = [&](double z, double y) { return uc(x[0], y); };
	function<double(double, double)> f2s = [&](double x, double z) { return us(x, y[0]); };
	function<double(double, double)> f2c = [&](double x, double z) { return uc(x, y[0]); };
	function<double(double, double)> f3s = [&](double z, double y) { return us(x[NX-1], y); };
	function<double(double, double)> f3c = [&](double z, double y) { return uc(x[NX-1], y); };
	function<double(double, double)> f4s = [&](double x, double z) { return us(x, y[NY-1]); };
	function<double(double, double)> f4c = [&](double x, double z) { return uc(x, y[NY-1]); };
	
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

	#pragma endregion
	#pragma region Задаем правую часть уравнения

	
	function<double(double, double)> Fs = [&](double x, double y) {return lambda * -1 * divgradS(x, y) - omega * sigma * uc(x, y) - pow(omega, 2) * hi * us(x, y); };
	function<double(double, double)> Fc = [&](double x, double y) {return lambda * -1 * divgradC(x, y) + omega * sigma * us(x, y) - pow(omega, 2) * hi * uc(x, y); };
	
	m.SetF(Fc, Fs);

	#pragma endregion


	#pragma region Вывод


#ifdef Print

	int pid = 12;

	m.Start();

	vector<double>	ResultU;
	vector<double>	ResultC;

	m.GetResult(ResultU, ResultC);

	for (size_t i = 0; i < (y.size() - 1)/2; i++)
		cout << setw(pid) << " ";
	cout << setw(pid) << "Us";
	for (size_t i = 0; i < y.size(); i++)
		cout << setw(pid) << " ";
	cout << setw(pid) << "Uc" << endl;

	for (size_t i = 0, i1 = 0, i2 = 0; i < y.size(); i++)
	{
		for (size_t j = 0; j < x.size(); j++, i1++)
		{
			cout << setw(pid) << ResultU[i1] << " ";
		}

		cout << "\t\t";

		for (size_t j = 0; j < x.size(); j++, i2++)
		{
			cout << setw(pid) << ResultC[i2] << " ";
		}
		cout << endl;
	}

	cout << endl;
	cout << "\t\t\tU" << endl;

	for (size_t i = 0, i1 = 0; i < y.size(); i++)
	{
		for (size_t j = 0; j < x.size(); j++, i1++)
		{
			cout << setw(pid) << ResultU[i1] * sin(0) + ResultC[i1] * cos(0);

		}
		cout << endl;
	}

	cout << endl;
	cout << "Погрешность" << endl;
	for (size_t i = 0; i < (y.size() - 1) / 2; i++)
		cout << setw(pid) << " ";
	cout << setw(pid) << "Us";
	for (size_t i = 0; i < y.size(); i++)
		cout << setw(pid) << " ";
	cout << setw(pid) << "Uc" << endl;

	double sum = 0;

	for (size_t i = 0, i1 = 0, i2 = 0; i < y.size(); i++)
	{
		for (size_t j = 0; j < x.size(); j++, i1++)
		{
			double p = us(x[j], y[i]) - ResultU[i1];
			cout << setw(pid) << p << " ";
			sum += pow(p, 2);
		}

		cout << "\t\t";

		for (size_t j = 0; j < x.size(); j++, i2++)
		{
			double p = uc(x[j], y[i]) - ResultC[i2];
			cout << setw(pid) << p << " ";
			sum += pow(p, 2);
		}
		cout << endl;
	}

	cout << "Норма вектора погрешности: " << sqrt(sum) << endl;
#endif // Print

#ifdef Check
	//ofstream out_t("mini2.txt");


	//out_t << setw(12) << "EPS" << setw(12) << "Labmda" << setw(12) << "Sigma" << setw(12) << "Omega" << setw(12) << "HI" << setw(12) << " MSG: " << setw(12) << "TIMELOS" << setw(12) << " IterLOS"<< setw(12) << " LU:" << setw(12) << "TIMELU" << endl;

	function<double(double, double)> Fs1;
	function<double(double, double)> Fc1;
	int pid = 12;

#ifdef CheckLambda
	{
	double eps = 1e-12, Sigma = 1e4, Omega = 1e2, HI = 1e-11;
	vector<double> Lambda = { 1e2, 1e3, 1e4, 1e5, 8e5 };
	double sum = 0;
	double LOSError = 0;
	double LUError = 0;
	double BSGError = 0;

	ofstream out_t("lambda.txt");
	out_t << "EPS LOS: " << eps << endl;
	out_t << "Sigma: " << Sigma << endl;
	out_t << "Omega: " << Omega << endl;
	out_t << "HI: " << HI << endl;
	out_t << setw(pid) << "Lambda " << setw(pid) << "LOS Time" << setw(pid) << "LOS Iter" << setw(pid) <<
		"LU Time" << setw(pid) << "BSG Time" << setw(pid) << "BSG Iter" << setw(pid) << "LOS error" << setw(pid) << "LU error" << setw(pid) << "BSG error"<< endl;
	for (size_t i = 0; i < Lambda.size(); i++)
	{
		Fs1 = [&](double x, double y) {return Lambda[i] * -1 * divgradS(x, y) - Omega * Sigma * uc(x, y) - pow(Omega, 2) * HI * us(x, y); };
		Fc1 = [&](double x, double y) {return Lambda[i] * -1 * divgradC(x, y) + Omega * Sigma * us(x, y) - pow(Omega, 2) * HI * uc(x, y); };

		m.SetF(Fc1, Fs1);
		m.SetLambda(Lambda[i]);
		m.SetSigma(Sigma);
		m.SetOmega(Omega);
		m.SetHi(HI);
		m.eps = eps;
						
		m.Start();

		sum = 0;

		for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
		{
			for (size_t j = 0; j < x.size(); j++, i1+=2)
				sum += pow(us(x[j], y[i]) - m.q1[i1], 2);
			for (size_t j = 0; j < x.size(); j++, i2+=2)
				sum += pow(uc(x[j], y[i]) - m.q1[i2], 2);
		}
		LOSError = sqrt(sum);
		sum = 0;

		for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
		{
			for (size_t j = 0; j < x.size(); j++, i1+=2)
				sum += pow(us(x[j], y[i]) - m.q2[i1], 2);

			for (size_t j = 0; j < x.size(); j++, i2+=2)
				sum += pow(uc(x[j], y[i]) - m.q2[i2], 2);
		}
		LUError = sqrt(sum);
		sum = 0;

		for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
		{
			for (size_t j = 0; j < x.size(); j++, i1 += 2)
				sum += pow(us(x[j], y[i]) - m.q[i1], 2);

			for (size_t j = 0; j < x.size(); j++, i2 += 2)
				sum += pow(uc(x[j], y[i]) - m.q[i2], 2);
		}

		BSGError = sqrt(sum);
		sum = 0;

		out_t << setw(pid) << Lambda[i] << setw(pid) << m.tq1 / 1000000.0 << setw(pid) << m.lastIter << setw(pid) << m.tq2 / 1000000.0 <<
			setw(pid) << m.tq3 / 1000000.0 << setw(pid) << m.lastIter2 << setw(pid) << LOSError << setw(pid) << LUError << setw(pid) << BSGError << endl;
	}
	out_t.close();
	}
#endif
#ifdef CheckSigma
	{
		double eps = 1e-12, Lambda = 1e4, Omega = 1e2, HI = 1e-11;
		vector<double> Sigma = { 0, 1, 1e2, 1e4, 1e6, 1e8 };
		double sum = 0;
		double LOSError = 0;
		double LUError = 0;
		double BSGError = 0;

		ofstream out_t("sigma.txt");
		out_t << "EPS LOS: " << eps << endl;
		out_t << "Lambda: " << Lambda << endl;
		out_t << "Omega: " << Omega << endl;
		out_t << "HI: " << HI << endl;
		out_t << setw(pid) << "Sigma " << setw(pid) << "LOS Time" << setw(pid) << "LOS Iter" << setw(pid) <<
			"LU Time" << setw(pid) << "BSG Time" << setw(pid) << "BSG Iter" << setw(pid) << "LOS error" << setw(pid) << "LU error" << setw(pid) << "BSG error" << endl;
		for (size_t i = 0; i < Sigma.size(); i++)
		{
			Fs1 = [&](double x, double y) {return Lambda * -1 * divgradS(x, y) - Omega * Sigma[i] * uc(x, y) - pow(Omega, 2) * HI * us(x, y); };
			Fc1 = [&](double x, double y) {return Lambda * -1 * divgradC(x, y) + Omega * Sigma[i] * us(x, y) - pow(Omega, 2) * HI * uc(x, y); };

			m.SetF(Fc1, Fs1);
			m.SetLambda(Lambda);
			m.SetSigma(Sigma[i]);
			m.SetOmega(Omega);
			m.SetHi(HI);
			m.eps = eps;

			m.Start();

			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q1[i1], 2);
				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q1[i2], 2);
			}
			LOSError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q2[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q2[i2], 2);
			}
			LUError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q[i2], 2);
			}

			BSGError = sqrt(sum);
			sum = 0;

			out_t << setw(pid) << Sigma[i] << setw(pid) << m.tq1 / 1000000.0 << setw(pid) << m.lastIter << setw(pid) << m.tq2 / 1000000.0 <<
				setw(pid) << m.tq3 / 1000000.0 << setw(pid) << m.lastIter2 << setw(pid) << LOSError << setw(pid) << LUError << setw(pid) << BSGError << endl;
		}
		out_t.close();
	}
#endif
#ifdef CheckOmega
	{
		double eps = 1e-12, Lambda = 1e4, Sigma = 1e4, HI = 1e-11;
		vector<double> Omega = { 1e-4, 1e-2, 1, 1e2, 1e4, 1e6, 1e8, 1e9 };
		double sum = 0;
		double LOSError = 0;
		double LUError = 0;
		double BSGError = 0;

		ofstream out_t("omega.txt");
		out_t << "EPS LOS: " << eps << endl;
		out_t << "Lambda: " << Lambda << endl;
		out_t << "Sigma: " << Sigma << endl;
		out_t << "HI: " << HI << endl;
		out_t << setw(pid) << "Omega " << setw(pid) << "LOS Time" << setw(pid) << "LOS Iter" << setw(pid) <<
			"LU Time" << setw(pid) << "BSG Time" << setw(pid) << "BSG Iter" << setw(pid) << "LOS error" << setw(pid) << "LU error" << setw(pid) << "BSG error" << endl;
		for (size_t i = 0; i < Omega.size(); i++)
		{
			Fs1 = [&](double x, double y) {return Lambda * -1 * divgradS(x, y) - Omega[i] * Sigma * uc(x, y) - pow(Omega[i], 2) * HI * us(x, y); };
			Fc1 = [&](double x, double y) {return Lambda * -1 * divgradC(x, y) + Omega[i] * Sigma * us(x, y) - pow(Omega[i], 2) * HI * uc(x, y); };

			m.SetF(Fc1, Fs1);
			m.SetLambda(Lambda);
			m.SetSigma(Sigma);
			m.SetOmega(Omega[i]);
			m.SetHi(HI);
			m.eps = eps;

			m.Start();

			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q1[i1], 2);
				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q1[i2], 2);
			}
			LOSError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q2[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q2[i2], 2);
			}
			LUError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q[i2], 2);
			}

			BSGError = sqrt(sum);
			sum = 0;

			out_t << setw(pid) << Omega[i] << setw(pid) << m.tq1 / 1000000.0 << setw(pid) << m.lastIter << setw(pid) << m.tq2 / 1000000.0 <<
				setw(pid) << m.tq3 / 1000000.0 << setw(pid) << m.lastIter2 << setw(pid) << LOSError << setw(pid) << LUError << setw(pid) << BSGError << endl;
		}
		out_t.close();
	}
#endif
#ifdef CheckHI
	{
		double eps = 1e-12, Lambda = 1e4, Sigma = 1e4, Omega = 1e2;
		vector<double> HI = { 8.81e-12, 1e-12, 1e-11, 1e-10 };
		double sum = 0;
		double LOSError = 0;
		double LUError = 0;
		double BSGError = 0;

		ofstream out_t("HI.txt");
		out_t << "EPS LOS: " << eps << endl;
		out_t << "Lambda: " << Lambda << endl;
		out_t << "Sigma: " << Sigma << endl;
		out_t << "Omega: " << Omega << endl;
		out_t << setw(pid) << "HI " << setw(pid) << "LOS Time" << setw(pid) << "LOS Iter" << setw(pid) <<
			"LU Time" << setw(pid) << "BSG Time" << setw(pid) << "BSG Iter" << setw(pid) << "LOS error" << setw(pid) << "LU error" << setw(pid) << "BSG error" << endl;
		for (size_t i = 0; i < HI.size(); i++)
		{
			Fs1 = [&](double x, double y) {return Lambda * -1 * divgradS(x, y) - Omega * Sigma * uc(x, y) - pow(Omega, 2) * HI[i] * us(x, y); };
			Fc1 = [&](double x, double y) {return Lambda * -1 * divgradC(x, y) + Omega * Sigma * us(x, y) - pow(Omega, 2) * HI[i] * uc(x, y); };

			m.SetF(Fc1, Fs1);
			m.SetLambda(Lambda);
			m.SetSigma(Sigma);
			m.SetOmega(Omega);
			m.SetHi(HI[i]);
			m.eps = eps;

			m.Start();

			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q1[i1], 2);
				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q1[i2], 2);
			}
			LOSError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q2[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q2[i2], 2);
			}
			LUError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q[i2], 2);
			}

			BSGError = sqrt(sum);
			sum = 0;

			out_t << setw(pid) << HI[i] << setw(pid) << m.tq1 / 1000000.0 << setw(pid) << m.lastIter << setw(pid) << m.tq2 / 1000000.0 <<
				setw(pid) << m.tq3 / 1000000.0 << setw(pid) << m.lastIter2 << setw(pid) << LOSError << setw(pid) << LUError << setw(pid) << BSGError << endl;
		}
		out_t.close();
	}
#endif

#ifdef CheckEPS
	{
		double Lambda = 1e4, Sigma = 1e4, Omega = 1e2, HI = 1e-11;
		vector<double> eps = { 1e-10, 1e-11, 1e-12, 1e-13, 1e-14 };
		double sum = 0;
		double LOSError = 0;
		double LUError = 0;
		double BSGError = 0;

		ofstream out_t("eps.txt");
		out_t << "Lambda: " << Lambda << endl;
		out_t << "Sigma: " << Sigma << endl;
		out_t << "Omega: " << Omega << endl;
		out_t << "HI: " << HI << endl;
		out_t << setw(pid) << "eps " << setw(pid) << "LOS Time" << setw(pid) << "LOS Iter" << setw(pid) <<
			"LU Time" << setw(pid) << "BSG Time" << setw(pid) << "BSG Iter" << setw(pid) << "LOS error" << setw(pid) << "LU error" << setw(pid) << "BSG error" << endl;
		for (size_t i = 0; i < eps.size(); i++)
		{
			Fs1 = [&](double x, double y) {return Lambda * -1 * divgradS(x, y) - Omega * Sigma * uc(x, y) - pow(Omega, 2) * HI * us(x, y); };
			Fc1 = [&](double x, double y) {return Lambda * -1 * divgradC(x, y) + Omega * Sigma * us(x, y) - pow(Omega, 2) * HI * uc(x, y); };

			m.SetF(Fc1, Fs1);
			m.SetLambda(Lambda);
			m.SetSigma(Sigma);
			m.SetOmega(Omega);
			m.SetHi(HI);
			m.eps = eps[i];

			m.Start();

			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q1[i1], 2);
				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q1[i2], 2);
			}
			LOSError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q2[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q2[i2], 2);
			}
			LUError = sqrt(sum);
			sum = 0;

			for (size_t i = 0, i1 = 0, i2 = 1; i < y.size(); i++)
			{
				for (size_t j = 0; j < x.size(); j++, i1 += 2)
					sum += pow(us(x[j], y[i]) - m.q[i1], 2);

				for (size_t j = 0; j < x.size(); j++, i2 += 2)
					sum += pow(uc(x[j], y[i]) - m.q[i2], 2);
			}

			BSGError = sqrt(sum);
			sum = 0;

			out_t << setw(pid) << eps[i] << setw(pid) << m.tq1 / 1000000.0 << setw(pid) << m.lastIter << setw(pid) << m.tq2 / 1000000.0 <<
				setw(pid) << m.tq3 / 1000000.0 << setw(pid) << m.lastIter2 << setw(pid) << LOSError << setw(pid) << LUError << setw(pid) << BSGError << endl;
		}
		out_t.close();
	}
#endif

#endif // Check



	#pragma endregion

	system("pause");
	return 0;
}