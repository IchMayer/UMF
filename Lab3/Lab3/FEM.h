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

		TypeBorder type;
		function<double(double, double)>& Us;
		function<double(double, double)>& Uc;

		//Задание первой краевой функции
		Border(function<double(double, double)> &fs, function<double(double, double)> &fc, TypeBorder t) : Us(fs), Uc(fc)
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

	void GetResult(vector<vector<double>> &Us, vector<vector<double>> &Uc)
	{
		Us = us;
		Uc = uc;
	}

	#pragma region SetFunctions
	//Задать базис
	void SetBasis(Basis b)
	{
		basis = b;
	}
	//Задать Сигму
	void SetSigma(double Sigma) { sigma = Sigma; }
	void SetOmega(double Omega) { omega = Omega; }
	void SetLambda(double l) { lambda = l; }
	//Добавление стороны
	void AddBorder(Border b)
	{
		border.push_back(b);
	}
	//Дабавление лямбды
	//Задание правой части
	void SetF(function<double(double, double, double)> func, function<double(double, double, double)> funs)
	{
		fc = func;
		fs = funs;
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
	//SetY
	void SetY(vector<double> Y)
	{
		y = Y;
	}
	#pragma endregion
private:
	#pragma region Privat struct
	struct ShapeIndex
	{
		int				ixMax;			//Максимальный индекс
		int				iyMax;
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
		void resize(int N)
		{
			n = N;
			x.resize(n);
			for (size_t i = 0; i < n; i++)
				x[i].resize(n);
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

	vector<double>	x;
	vector<double>	y;

	double			lambda;			//Лямбда
	double			sigma;			//Сигма
	double			omega;			//Омега
	function<double(double, double, double)> fc;
	function<double(double, double, double)> fs;

	//Результат
	vector<vector<double>>	us;		//Искомая функция
	vector<vector<double>>	uc;		//Искомая функция

	//Переходные вычисления
	vector<vector<double>> A;		//Матрица левой части в диагональном виде
	vector<double>	b;				//Вектора правой части
	vector<double>	q;			

	double _x, _xl;
	double _y, _yl;

	Matrix P;
	Matrix C;

	Matrix M;
	Matrix G;
	Matrix cF;

	ShapeIndex		index;

	//Обработка введеных данных
	void Preprocessing()
	{
		index.ixMax = x.size();
		index.iyMax = y.size();

		int m = index.ixMax * index.iyMax;
		A.resize(m);
		for (size_t i = 0; i < m; i++)
			A[i].resize(m);

		b.resize(m);

		int order = basis.order + 1;

		P.resize(order);
		C.resize(order);

		if (!q.size())
			q.resize(m);

		//k = (index.ixMax - 1.0) / 2.0;	//Определяем количество локальных матриц
		
		CreateMG();
	}
	void CreateMG()
	{
		switch (basis.type)
		{
		case Basis::Lagrange:
			switch (basis.order)
			{
			case 1:
				G.x = vector<vector<double>>{ {2, 1, -2, -1}, {1, 2, -1, -2}, {-2, -1, 2, 1}, {-1, -2, 1, 2} };
				G = G *(this->lambda / 6.0);
				G.n = 4;
				M.x = vector<vector<double>>{ {4, 2, 2, 1}, {2, 4, 1, 2}, {2, 1, 4, 2}, {1, 2, 2, 4}};
				cF.x = M.x;
				M = M * (this->sigma / 36.0);
				M.n = 4;
				cF = cF/ 36.0;
				cF.n = 4;
				break;
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
			for (size_t i = 0; i < q.size(); i++)
			{
				for (size_t j = 0; j < q.size(); j++)
					A[i*index.ixMax][j] = 0;
				A[i * index.ixMax][i * index.ixMax] = 1;

				b[i * index.ixMax] = border[0].U(_x, _y);
			}

			break;	
		default:
			break;
		}

		switch (border[1].type)
		{
		case Border::First:
			for (size_t i = 0; i < q.size(); i++)
			{
				for (size_t j = 0; j < q.size(); j++)
					A[i][j] = 0;
				A[i][i] = 1;

				b[i] = border[1].U(_x, _y);
			}

			break;
		default:
			break;
		}

		switch (border[2].type)
		{
		case Border::First:
			for (size_t i = 0; i < q.size(); i++)
			{
				for (size_t j = 0; j < q.size(); j++)
					A[(i + 1)* index.ixMax - 1][j] = 0;
				A[(i + 1)* index.ixMax - 1][(i + 1)* index.ixMax - 1] = 1;

				b[(i + 1)* index.ixMax - 1] = border[2].U(_x, _y);
			}

			break;
		default:
			break;
		}

		switch (border[3].type)
		{
		case Border::First:
			for (size_t i = 0; i < q.size(); i++)
			{
				for (size_t j = 0; j < q.size(); j++)
					A[(i + 1)* index.ixMax - 1][j] = 0;
				A[(i + 1)* index.ixMax - 1][(i + 1)* index.ixMax - 1] = 1;

				b[(i + 1)* index.ixMax - 1] = border[1].U(_x, _y);
			}

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
