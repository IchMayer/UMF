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

	//Форма
	struct Shape
	{
		enum ShapeType
		{
			Recrangle,
			T,
			Line
		};

		Shape(int N = 2, double xmax = 1, ShapeType t = Recrangle ,double xmin = 0)
		{
			n = N;
			type = t;
			xMin = xmin;
			xMax = xmax;
		}
			
		ShapeType type;
		double xMin, xMax;
		int n;
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

//#pragma region 2d
//	//Граничные условия
//	struct Border
//	{
//	public:
//		enum TypeBorder
//		{
//			First,
//			Second,
//			Third
//		};
//
//		struct Grad
//		{
//			Grad(function<double(double, double)> a, function<double(double, double)> b) :x(a), y(b) {}
//			function<double(double, double)>& x;
//			function<double(double, double)>& y;
//		};
//
//		TypeBorder type;
//		function<double(double, double)>& U;
//		Grad Thetta;
//
//		//Задание первой краевой функции
//		Border(function<double(double, double)> &f) : U(f), Thetta(f, f)
//		{
//			type = First;
//		}
//		//Задание второй или третьей краевой функции
//		Border(function<double(double, double)> &fx, function<double(double, double)> &fy, TypeBorder type = Second) :Thetta(fx, fy), U(fx)
//		{
//			this->type = type;
//		}
//	};
//#pragma endregion

	FEM() { lambda = 1; };
	~FEM() {};

	double eps;
	int maxIter;
	int lastIter;

	void Start()
	{
		checkdata();

		Preprocessing();

		auto hinstLib = LoadLibrary(TEXT("LDU.dll"));
		Gause LDU = (Gause)GetProcAddress(hinstLib, "Gause");

		vector<double> localF(basis.order + 1);
		auto qlast = q;

		for (size_t i = 0; i < nt; i++)
		{
			t = t0 + i * ht;
			lastIter = 0;
			do
			{
				ClearLastIter();
				for (size_t i = 0; i < k; i++)
				{
					for (size_t j = 0; j <= basis.order; j++)
						localF[j] = f(shape.xMin + hx * i, t);
					blocal = (CreateMatrixС(hx) + M) * localF;

					//!!!! сдлеать сигму
					M = CreateMatrixС(hx * gamma) * ht;		//Равномераная сетка по времени
					Alocal = M + G;

					AddLocalMatrix(i);
				}

				AddBorderInMatrix();

				if (Norm((A * q) - b) / Norm(b) < eps)	break;

				qlast = q;
				//Тут решаем уравнение A q = b
				LDU(A, b, q);

			} while (++lastIter < maxIter && Norm(q - qlast) / Norm(q) > eps);

			//т.к. базис лагранджа то u = q
			u.push_back(q);
		}
	}

	void Clear()
	{

	}

	//Задать базис
	void SetBasis(Basis b)
	{
		basis = b;
	}

	//Задать форму
	void SetShape(Shape s)
	{
		shape = s;
	}

	//Задать начатьное приближение
	void SetQ0(vector<double> q0)
	{
		q = q0;
	}

	//Задать Сигму
	void SetSigma(function<double(double)> Sigma)
	{
		sigma = Sigma;
	}

	//Задать Гамму
	void SetGamma(double Gamma)
	{
		gamma = Gamma;
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

	//Set t0
	void SetT0(double t)
	{
		t0 = t;
	}
	//Set tn
	void SetNt(double t)
	{
		nt = t;
	}
	//Set ht
	void SetHt(double t)
	{
		ht = t;
	}

private:
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
		//Matrix operator+(Matrix b)
		//{
		//	return Matrix(vector<vector<double>>{ {7, -8, 1}, { -8, 16, -5 }, { 1, -8, 7 }});
		//}
	};

	//Дано
	Basis			basis;			//Базис		
	Shape			shape;			//Форма в общих коорлинатах
	vector<Border>	border;			//Границы
	vector<double>	q;				//Начальное значение
	double			t0;				//Начальное время
	double			nt;				//Количество шагов
	double			ht;				//Шаг по времени
	double			lambda;			//Лямбда
	double			gamma;			//Гамма
	function<double(double)> sigma;	//Сигма
	function<double(double, double)> f;		//Правая часть 

	//Результат
	vector<vector<double>>	u;		//Искомая функция

	//Переходные вычисления
	double			hx;				//Шаг по x
	ShapeIndex		index;			//Форма в координатах интексов
	int				k;				//Количество локальных матриц
	double			t;				//Текущее время

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
		hx = (shape.xMax - shape.xMin) / shape.n;
		index.ixMax = shape.n;

		A.resize(shape.n);
		for (size_t i = 0; i < shape.n; i++)
			A[i].resize(shape.n);

		b.resize(shape.n);
		blocal.resize(basis.order + 1);
		Alocal = Matrix(basis.order + 1);

		if (!q.size())
			q.resize(shape.n);

		k = (shape.n - 1.0) / 2.0;	//Определяем количество локальных матриц
		
		add_basis(basis);
	}

	//Добавление базисныйх функций 
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


				double c = lambda / 3 / 2 / hx;
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
		if (ht == 0.0)
			throw "error ht";
		if (lambda == 0)
			throw "error lambda";
	}

	//Создание матрицы массы
	Matrix CreateMatrixС(vector<double> vx = {1, 1, 1})
	{
		Matrix C(3);
		C.x[0][0] = vx[0] * 0.09285714285714286 +	vx[1] * 0.04761904761904761 +	vx[2] * -7.1428571428571415e-3;
		C.x[0][1] = vx[0] * 0.04761904761904761 +	vx[1] * 0.0380952380952381 +	vx[2] * -0.01904761904761905;
		C.x[0][2] = vx[0] * -7.1428571428571415e-3 + vx[1]* -0.01904761904761905 +	vx[2] * -7.1428571428571415e-3;
		C.x[1][0] = vx[0] * 0.04761904761904761 +	vx[1] * 0.0380952380952381 +	vx[2] * -0.01904761904761905;
		C.x[1][1] = vx[0] * 0.0380952380952381 +	vx[1] * 0.4571428571428571 +	vx[2] * 0.03809523809523809;
		C.x[1][2] = vx[0] * -0.01904761904761905 +	vx[1] * 0.03809523809523809 +	vx[2] * 0.04761904761904761;
		C.x[2][0] = vx[0] * -7.1428571428571415e-3 + vx[1] * -0.01904761904761905 + vx[2] * -7.1428571428571415e-3;
		C.x[2][1] = vx[0] * -0.01904761904761905 +	vx[1] * 0.03809523809523809 +	vx[2] * 0.04761904761904761;
		C.x[2][2] = vx[0] * -7.1428571428571415e-3 + vx[1] * 0.04761904761904761 +	vx[2] * 0.09285714285714286;
		return C;
	}

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

	void AddBorderInMatrix()
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

	double Norm(vector<double> b)
	{
		double sum = 0;
		for (auto &i : b)
			sum = i * pow(i, 2);
		return sqrt(sum);
	}
};