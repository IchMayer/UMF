#include "Includes.h"

typedef void(*Gause)(std::vector<std::vector<double>> &A, std::vector<double> &F, std::vector<double> &X);

vector<double> operator*(vector<vector<double>> a, vector<double> b)
{
	vector<double> c(b.size());
	double sum = 0;
	int m = a.size();
	for (size_t i = 0; i < m; i++)
	{
		sum = 0;
		for (size_t j = 0; j < m; j++)
			sum += a[i][j] * b[j];
		c[i] = sum;
	}
	return c;
}
vector<double> operator+(vector<double> a, vector<double> b)
{
	vector<double> c(a.size());
	double sum = 0;
	for (size_t i = 0; i < a.size(); i++)
		c[i] = a[i] + b[i];
	return c;
}
vector<double> operator-(vector<double> a, vector<double> b)
{
	vector<double> c(a.size());
	double sum = 0;
	for (size_t i = 0; i < a.size(); i++)
		c[i] = a[i] - b[i];
	return c;
}
vector<double> operator*(vector<double> a, double b)
{
	for (size_t i = 0; i < a.size(); i++)
		a[i] *= b;
	return a;
}

class LUsolve
{
public:
	LUsolve() {};
	~LUsolve() {};
	
	vector<vector<double>> A;
	vector<double> F;
	vector<double> X;
	int n;	//Размерность матрицы А
	int m;	//Количество лент
	int d;	//Номер главной диаганали

private:

	vector<vector<double>> LU;

	void CountLU()
	{
		m = A.size();
		n = A[0].size();

		LU.resize(A.size());
		for (size_t i = 0; i < m; i++)
			LU[i].resize(n);

		double sum;

		for (size_t i = 0; i < n; i++)
		{

			for (size_t j = max(0, d - i); j < d; j++)
			{
				sum = 0;
				for (size_t k = 0; k < i; k++)
					sum += LU[i][k] * LU[k][j];
				LU[j][i] = A[j][i] - sum;
			}
		}
	}
	void CountX()
	{

	}
};

class FEM
{
public:
	#pragma region Public struct
	//Базисы
	struct Basis
	{
		enum BasisType
		{
			Lagrange,
			Hermite
		};

		BasisType type;
		int order;
	};
	//Граничные условия
	struct Border
	{
	public:
		enum TypeBorder
		{
			First,
			Second,
			Third
		};

		struct Grad
		{
			Grad(function<double(double)> a) :x(a){}
			function<double(double)>& x;
		};

		TypeBorder type;
		function<double(double)>& U;
		Grad Thetta;

		//Задание первой краевой функции
		Border(function<double(double)> &f, TypeBorder t) : U(f), Thetta(f)
		{
			type = t;
		}
	};
	#pragma endregion

	FEM() { lambda = 1;};
	~FEM() {};

	double eps;			//Погрешнасть решения
	int maxIter;		//Максимальное еоличестов итераций на один слой
	int lastIter;		//Текущее количество итераций на слое

	vector<int>		It;

	//Запуск метода конечных элементов
	void Start()
	{
		checkdata();

		Preprocessing();

		auto hinstLib = LoadLibrary(TEXT("LDU.dll"));
		Gause LDU = (Gause)GetProcAddress(hinstLib, "Gause");

		vector<double> Zsigma(basis.order + 1);
		vector<double> Zsigmalast(basis.order + 1);
		vector<double> localF(basis.order + 1);
		vector<double> qlast = q;
		vector<double> qlast2 = q;
		double len;
		double dt;
		Matrix c;
		Matrix Alast;

		for (size_t l = 1; l < t.size(); l++)
		{
			lastIter = 0;
			do
			{
				dt = t[l] - t[l - 1];
				ClearLastIter();
				for (size_t i = 0; i < x.size() - 1; i++)
				{
					localF[0] = f(x[i], t[l]);
					localF[1] = f((x[i + 1] + x[i]) / 2, t[l]);
					localF[2] = f(x[i + 1], t[l]);

					len = x[i + 1] - x[i];

					Zsigma[0] = sigma(x[i], t[l], (q[2 * i + 1] - q[2 * i]) / (len / 2));
					Zsigma[1] = sigma((x[i+1]+x[i])/2, t[l], (q[2 * i + 2] - q[2 * i]) / len);
					Zsigma[2] = sigma(x[i+1], t[l], (q[2 * i + 2] - q[2 * i + 1]) / (len / 2));

					M = CreateMatrixС(Zsigma * len) / dt;

					blocal = (CreateMatrixС(len)) * localF + M * vector<double>{ qlast2[2 * i], qlast2[2 * i + 1], qlast2[2 * i + 2] };

					Alocal = M + G / len;

					AddLocalMatrix(i);
				}

				AddBorderInMatrix(t[l]);

				auto kj = Norm((A * q) - b);

				if (Norm((A * q) - b) / Norm(b) < eps)	break;

				qlast = q;
				//Тут решаем уравнение A q = b
				LDU(A, b, q);

			} while (++lastIter < maxIter && Norm(q - qlast) / Norm(q) > eps);

			It.push_back(lastIter);

			qlast2 = q;
			//т.к. базис лагранджа то u = q
			u.push_back(q);
		}
	}

	//Высвобождение памяти //кототрого нет!!!
	void Clear()
	{

	}

	void GetResult(vector<vector<double>> &U)
	{
		U = u;
	}

	#pragma region SetFunctions
	//Задать базис
	void SetBasis(Basis b)
	{
		basis = b;
	}
	//Задать начатьное приближение
	void SetQ0(vector<double> q0)
	{
		q = q0;
	}
	//Задать Сигму
	void SetSigma(function<double(double, double, double)> Sigma)
	{
		sigma = Sigma;
	}
	//Добавление стороны
	void AddBorder(Border b)
	{
		border.push_back(b);
	}
	//Дабавление лямбды
	void SetLambda(double l)
	{
		lambda = l;
	}
	//Задание правой части
	void SetF(function<double(double, double)> func)
	{
		f = func;
	}
	//Set t
	void SetT(vector<double> T)
	{
		t = T;
	}
	//SetX
	void SetX(vector<double> X)
	{
		x = X;
	}
	#pragma endregion
private:
	#pragma region Privat struct
	struct ShapeIndex
	{
		int				ixMax;			//Максимальный индекс			
	};
	
	struct Matrix
	{
		Matrix()
		{
			n = 0;
		}
		Matrix(int N)
		{
			n = N;
			x.resize(n);
			for (size_t i = 0; i < n; i++)
				x[i].resize(n);
		}
		Matrix(vector<vector<double>> x0)
		{
			x = x0;
			n = x0.size();
		}

		vector<vector<double>> x;
		size_t n;

		void Clear()
		{
			x.clear();
		}

		Matrix operator*(double b)
		{
			Matrix a(x);
			for (size_t i = 0; i < a.x.size(); i++)
				for (size_t j = 0; j < a.x[i].size(); j++)
					a.x[i][j] *= b;
			return a;
		}
		Matrix operator/(double b)
		{
			Matrix a(x);
			for (size_t i = 0; i < a.x.size(); i++)
				for (size_t j = 0; j < a.x[i].size(); j++)
					a.x[i][j] /= b;
			return a;
		}
		Matrix operator+(Matrix b)
		{
			Matrix c(3);
			for (size_t i = 0; i < n; i++)
				for (size_t j = 0; j < n; j++)
					c.x[i][j] = x[i][j] + b.x[i][j];
			return c;
		}
		void operator*=(double b)
		{
			for (size_t i = 0; i < x.size(); i++)
				for (size_t j = 0; j < x[i].size(); j++)
					x[i][j] *= b;
		}
		vector<double> operator*(vector<double> b)
		{
			vector<double> a(n);

			for (size_t i = 0; i < n; i++)
				for (size_t j = 0; j < n; j++)
					a[i] += x[i][j] * b[j];

			return a;
		}
	};
	#pragma endregion

	//Дано
	Basis			basis;			//Базис		
	vector<Border>	border;			//Границы

	vector<double>  t;				//Время
	vector<double>	x;				//Координаты точке
	vector<double>	q;				//Начальное значение

	double			lambda;			//Лямбда
	function<double(double, double, double)> sigma;	//Сигма
	function<double(double, double)> f;		//Правая часть 

	//Результат
	vector<vector<double>>	u;		//Искомая функция

	//Переходные вычисления
	ShapeIndex		index;			//Форма в координатах интексов
	int				k;				//Количество локальных матриц

	Matrix			M;				//Матрица массы
	Matrix			G;				//Матрица жесткости
	Matrix			Alocal;			//Локальная матрица левой части

	//vector<function<double(double)>> fbasis;		//Базисы
	vector<vector<double>> A;		//Матрица левой части в диагональном виде
	vector<double>	b;				//Вектора правой части
	vector<double>	blocal;			//Локальный вектор правой части


	//Обработка введеных данных
	void Preprocessing()
	{
		index.ixMax = 2 * x.size() - 1;

		A.resize(index.ixMax);
		for (size_t i = 0; i < index.ixMax; i++)
			A[i].resize(index.ixMax);

		b.resize(index.ixMax);
		blocal.resize(basis.order + 1);
		Alocal = Matrix(basis.order + 1);

		if (!q.size())
			q.resize(index.ixMax);

		k = (index.ixMax - 1.0) / 2.0;	//Определяем количество локальных матриц
		
		add_basis(basis);
	}
	//Добавление базисныйх функций и матрыцы жесткости
	void add_basis(Basis b)
	{
		switch (b.type)
		{
		case b.Lagrange:
			switch (b.order)
			{
			case 2:
			{
				//fbasis.resize(3);
				//fbasis[0] = [](double x) { return 2 * (x - 0.5) * (x - 1); };
				//fbasis[1] = [](double x) { return -4 * x  * (x - 1); };
				//fbasis[2] = [](double x) { return 2 * x * (x - 0.5); };


				double c = lambda / 3;
				G = Matrix({ {7, -8, 1}, {-8, 16, -8}, {1, -8, 7} });
				G *= c;

				M.x.resize(3);
				for (size_t i = 0; i < M.x.size(); i++)
					M.x[i].resize(3);
				break;
			}
			default:
				break;
			}
			break;
		default:
			break;
		}
	}
	//Провкрка введенных данных
	void checkdata()
	{
		if (border.size() < 2)
			throw "all boundaries are not indicated";
		if (lambda == 0)
			throw "error lambda";
	}
	//Создание матрицы массы
	Matrix CreateMatrixС(vector<double> vx)
	{
		Matrix C(3);
		C.x[0][0] = vx[0] * 0.09285714285714286		+	vx[1] * 0.04761904761904761		+	vx[2] * -7.1428571428571415e-3;
		C.x[0][1] = vx[0] * 0.04761904761904761		+	vx[1] * 0.0380952380952381		+	vx[2] * -0.01904761904761905;
		C.x[0][2] = vx[0] * -7.1428571428571415e-3	+	vx[1]* -0.01904761904761905		+	vx[2] * -7.1428571428571415e-3;
		C.x[1][0] = vx[0] * 0.04761904761904761		+	vx[1] * 0.0380952380952381		+	vx[2] * -0.01904761904761905;
		C.x[1][1] = vx[0] * 0.0380952380952381		+	vx[1] * 0.4571428571428571		+	vx[2] * 0.03809523809523809;
		C.x[1][2] = vx[0] * -0.01904761904761905	+	vx[1] * 0.03809523809523809		+	vx[2] * 0.04761904761904761;
		C.x[2][0] = vx[0] * -7.1428571428571415e-3	+	vx[1] * -0.01904761904761905	+	vx[2] * -7.1428571428571415e-3;
		C.x[2][1] = vx[0] * -0.01904761904761905	+	vx[1] * 0.03809523809523809		+	vx[2] * 0.04761904761904761;
		C.x[2][2] = vx[0] * -7.1428571428571415e-3	+	vx[1] * 0.04761904761904761		+	vx[2] * 0.09285714285714286;

		return C;
	}
	//Создание матрицы массы
	Matrix CreateMatrixС(double x)
	{
		Matrix C(3);
		C.x[0][0] = x * 0.13333333333333333;
		C.x[0][1] = x * 0.06666666666666667;
		C.x[0][2] = x * -0.03333333333333333;
		C.x[1][0] = x * 0.06666666666666667;
		C.x[1][1] = x * 0.5333333333333333;
		C.x[1][2] = x * 0.06666666666666667;
		C.x[2][0] = x * -0.03333333333333333;
		C.x[2][1] = x * 0.06666666666666667;
		C.x[2][2] = x * 0.13333333333333333;
		return C;
	}
	//Добавление локальных матриц в глобальную
	void AddLocalMatrix(int location)
	{
		for (size_t i = 0, il = location * basis.order; i <= basis.order; i++, il++)
		{
			for (size_t j = 0, jl = location * basis.order; j <= basis.order; j++, jl++)
				A[il][jl] += Alocal.x[i][j];
			b[il] += blocal[i];
		}
	}
	//Очистка A и b c предыдущей итерации
	void ClearLastIter()
	{
		for (size_t i = 0; i < A.size(); i++)
		{
			for (size_t j = 0; j < A[i].size(); j++)
				A[i][j] = 0;
			b[i] = 0;
		}
	}
	//Добавление краевых условаий в матрицу А и b
	void AddBorderInMatrix(double t)
	{

		switch (border[0].type)
		{
		case Border::First:
			A[0][0] = 1;
			for (size_t i = 1; i < A[0].size(); i++)
				A[0][i] = 0;
			b[0] = border[0].U(t);

			break;	
		default:
			break;
		}

		switch (border[1].type)
		{
		case Border::First:
			A[A.size() - 1][A.size() - 1] = 1;
			for (size_t i = 0; i < A[0].size() - 1; i++)
				A[A.size() - 1][i] = 0;
			b[A.size() - 1] = border[1].U(t);

			break;
		default:
			break;
		}
	}
	//Норма
	double Norm(vector<double> b)
	{
		double sum = 0;
		for (auto &i : b)
			sum += pow(i, 2);
		return sqrt(sum);
	}
};
