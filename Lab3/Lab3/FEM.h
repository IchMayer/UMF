#include "Includes.h"

typedef void(*Gause)(std::vector<std::vector<double>> &A, std::vector<double> &F, std::vector<double> &X);
typedef size_t(*MSGSolver)(std::vector<size_t> &iig, std::vector<size_t> &ijg, std::vector<double> &gl, std::vector<double> &gu, std::vector<double> &di, std::vector<double> &F, std::vector<double> &X, double EPS);
typedef void(*LDU)(std::vector<double> al, std::vector<double> au, std::vector<double> di, std::vector<size_t> &ia, std::vector<double> &F, std::vector<double> &X);

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
	int lastIter2;		//Текущее количество итераций на слое

	long long tq1, tq2, tq3;
	vector<double>	q1;
	vector<double>	q2;
	vector<double>	q;

	//Запуск метода конечных элементов
	void Start()
	{
		checkdata();

		Preprocessing();

		auto hinstLib = LoadLibrary(TEXT("LDU.dll"));
		Gause gause = (Gause)GetProcAddress(hinstLib, "Gause");
		LDU ldu = (LDU)GetProcAddress(hinstLib, "LDU");
		MSGSolver MSG_LOS_BCGSTAB = (MSGSolver)GetProcAddress(hinstLib, "MSGSolver");
		MSGSolver LOS = (MSGSolver)GetProcAddress(hinstLib, "LOS");
		MSGSolver BSG_STAB = (MSGSolver)GetProcAddress(hinstLib, "BSG_STAB");

		vector<double> localFs(pow(basis.order + 1, 2));
		vector<double> localFc(pow(basis.order + 1, 2));


		double hx, hy;
		vector<double> localbs(pow(basis.order + 1, 2));
		vector<double> localbc(pow(basis.order + 1, 2));

		lastIter = 0;
		ClearLastIter();
		for (size_t i = 0; i < x.size() - 1; i++)
		{
			for (size_t j = 0; j < y.size() - 1; j++)
			{
				localFs[0] = fs(x[i], y[j]);
				localFs[1] = fs(x[i + 1], y[j]);
				localFs[2] = fs(x[i], y[j + 1]);
				localFs[3] = fs(x[i + 1], y[j + 1]);

				localFc[0] = fc(x[i], y[j]);
				localFc[1] = fc(x[i + 1], y[j]);
				localFc[2] = fc(x[i], y[j + 1]);
				localFc[3] = fc(x[i + 1], y[j + 1]);

				hx = x[i + 1] - x[i];
				hy = y[j + 1] - y[j];

				P = Gx * (hy / hx) + Gy * (hx / hy) - M * (pow(omega, 2) * hi * hx * hy);
				C = M * (omega * sigma * hx * hy);

				CreateLocalA();
					
				localbc = (M *(hx * hy)) * localFc;
				localbs = (M *(hx * hy)) * localFs;

				for (size_t i = 0; i < localFc.size(); i++)
				{
					blocal[2 * i] = localbs[i];
					blocal[2 * i + 1] = localbc[i];
				}
				
				AddLocalMatrix(j * index.ixMax + i, (j+1)*index.ixMax + i);
			}
		}

		AddBorderInMatrix();
		
		q1.resize(q.size());
		q2.resize(q.size());

		Prof2Raz();
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		//lastIter = MSG(ig, jg, ggl, ggu, di, b, q1, eps);
		lastIter = LOS(ig, jg, ggl, ggu, di, b, q1, eps);
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		tq1 = duration_cast<microseconds>(t2 - t1).count();

		di.clear();
		ggu.clear();
		ggl.clear();
		ig.clear();
		jg.clear();

		FULL2Prof();
		t1 = high_resolution_clock::now();
		ldu(ggl, ggu, di, ig, b, q2);
		t2 = high_resolution_clock::now();
		tq2 = duration_cast<microseconds>(t2 - t1).count();

		di.clear();
		ggu.clear();
		ggl.clear();
		ig.clear();
		jg.clear();

		Prof2Raz();
		t1 = high_resolution_clock::now();
		//lastIter2 = BSG_STAB(ig, jg, ggl, ggu, di, b, q, eps);
		t2 = high_resolution_clock::now();
		tq3 = duration_cast<microseconds>(t2 - t1).count();

		di.clear();
		ggu.clear();
		ggl.clear();
		ig.clear();
		jg.clear();

		gause(A, b, q);

#ifdef _DEBUG
		system("cls");
		for (size_t i = 0; i < A.size(); i++)
		{
			for (size_t j = 0; j < A[i].size(); j++)
			{
				cout << setw(10) << A[i][j] << " ";
			}
			cout << "\t\t\t" << b[i] << endl;
		}

		for (size_t i = 0; i < q.size(); i++)
		{
			cout << q[i] << "   "<< q1[i] << "   "<< q2[i] << endl;
		}

#endif // _DEBUG

		//т.к. базис лагранджа то u = q
		for (size_t i = 0; i < q.size() / 2; i++)
		{
			us[i] = q[2 *i];
			uc[i] = q[2 *i + 1];
		}
	}

	//Высвобождение памяти //кототрого нет!!!
	void Clear()
	{

	}

	void GetResult(vector<double> &Us, vector<double> &Uc)
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
	void SetHi(double Hi) { hi = Hi; }
	//Добавление стороны
	void AddBorder(Border b)
	{
		border.push_back(b);
	}
	//Дабавление лямбды
	//Задание правой части
	void SetF(function<double(double, double)> func, function<double(double, double)> funs)
	{
		fc = func;
		fs = funs;
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
			Matrix c(x);
			for (size_t i = 0; i < n; i++)
				for (size_t j = 0; j < n; j++)
					c.x[i][j] += b.x[i][j];
			return c;
		}
		Matrix operator-(Matrix b)
		{
			Matrix c(x);
			for (size_t i = 0; i < n; i++)
				for (size_t j = 0; j < n; j++)
					c.x[i][j] -= b.x[i][j];
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
	double			hi;
	function<double(double, double)> fc;
	function<double(double, double)> fs;

	//Результат
	vector<double>	us;		//Искомая функция
	vector<double>	uc;		//Искомая функция

	//Переходные вычисления
	int m;
	vector<vector<double>> A;		//Матрица левой части в диагональном виде
	vector<double>	b;				//Вектора правой части				

	Matrix Alocal;
	vector<double> blocal;

	Matrix P;
	Matrix C;

	Matrix M;
	Matrix Gx;			//hy/hx
	Matrix Gy;			//hx/hy

	vector<double> di;	//Диагональные элементы матрицы А
	vector<double> ggu;	//Верхний треугольник матрицы А в разреженном формате
	vector<double> ggl;	//Нижний треугольний матрицы А в разреженном формате
	vector<size_t> ig;// Массив индексов
	vector<size_t> jg; // Другой массив индексов

	ShapeIndex		index;

	void Prof2Raz()
	{
		ig.push_back(0);

		for (size_t i = 0, j1 = 0; i < m; i++)
		{
			for (size_t j = 0; j < i; j++)
			{
				if (A[i][j] || A[j][i])
				{
					ggu.push_back(A[j][i]);
					ggl.push_back(A[i][j]);
					jg.push_back(j);
					j1++;
				}
			}
			ig.push_back(j1);
			di.push_back(A[i][i]);
		}

	}

	void FULL2Prof()
	{
		di.clear();
		ggu.clear();
		ggl.clear();
		ig.clear();
		jg.clear();
		ig.push_back(0);

		for (size_t i = 0, j1 = 0; i < m; i++)
		{
			for (size_t j = 0, b = 0; j < i; j++)
			{
				if (A[i][j] || A[j][i] || b)
				{
					ggu.push_back(A[j][i]);
					ggl.push_back(A[i][j]);
					j1++;
					b = 1;
				}
			}
			ig.push_back(j1);
			di.push_back(A[i][i]);
		}

	}

	//Обработка введеных данных
	void Preprocessing()
	{
		index.ixMax = x.size();
		index.iyMax = y.size();

		m = 2 * index.ixMax * index.iyMax;
		
		A.resize(m);
		for (size_t i = 0; i < m; i++)
			A[i].resize(m);

		b.resize(m);
		q.resize(m);


		int order = basis.order + 1;

		Alocal.resize(2 * pow(order, 2));
		blocal.resize(2 * pow(order, 2));

		P.resize(pow(order, 2));
		C.resize(pow(order, 2));
		
		CreateMG();
		us.resize(m / 2);
		uc.resize(m / 2);
	}
	void CreateMG()
	{
		Gy.x = vector<vector<double>>{ {2, 1, -2, -1}, {1, 2, -1, -2}, {-2, -1, 2, 1}, {-1, -2, 1, 2} };
		Gy = Gy *(this->lambda / 6.0);
		Gy.n = 4;

		Gx.x = vector<vector<double>>{ {2, -2, 1, -1}, {-2, 2, -1, 1}, {1, -1, 2, -2}, {-1, 1, -2, 2} };
		Gx = Gx * (this->lambda / 6.0);
		Gx.n = 4;

		M.x = vector<vector<double>>{ {4, 2, 2, 1}, {2, 4, 1, 2}, {2, 1, 4, 2}, {1, 2, 2, 4}};
		M = M / 36.0;
		M.n = 4;
	}
	void CreateLocalA()
	{
		for (size_t i = 0; i < P.n; i++)
			for (size_t j = 0; j < P.n; j++)
				Alocal.x[2 * i + 1][2 * j + 1] = Alocal.x[2 * i][2 * j] = P.x[i][j];

		for (size_t i = 0; i < C.n; i++)
			for (size_t j = 0; j < C.n; j++)
			{
				Alocal.x[2 * i][2 * j + 1] = -C.x[i][j];
				Alocal.x[2 * i + 1][2 * j] = C.x[i][j];
			}
	}
	//Провкрка введенных данных
	void checkdata()
	{

	}
	//Добавление локальных матриц в глобальную
	void AddLocalMatrix(int x0, int x1)
	{
		for (size_t i = 0, il = 2 * x0 * basis.order; i < 2 * (basis.order + 1); i++, il++)
		{
			for (size_t j = 0, jl = 2 * x0 * basis.order; j < 2 * (basis.order + 1); j++, jl++)
				A[il][jl] += Alocal.x[i][j];
		}

		for (size_t i = 0, il = 2 * x0 * basis.order; i < 2 * (basis.order + 1); i++, il++)
		{
			for (size_t j = 0, j2 = 2 * (basis.order + 1), jl = 2 * x1 * basis.order; j < 2 * (basis.order + 1); j++, jl++, j2++)
				A[il][jl] += Alocal.x[i][j2];
		}

		for (size_t i = 0, i1 = 2 * (basis.order + 1), il = 2 * x1 * basis.order; i < 2 * (basis.order + 1); i++, il++, i1++)
		{
			for (size_t j = 0, jl = 2 * x0 * basis.order; j < 2 * (basis.order + 1); j++, jl++)
				A[il][jl] += Alocal.x[i1][j];
		}

		for (size_t i = 0, i1 = 2 * (basis.order + 1), il = 2 * x1 * basis.order; i < 2 * (basis.order + 1); i++, il++, i1++)
		{
			for (size_t j = 0, j1 = 2 * (basis.order + 1), jl = 2 * x1 * basis.order; j < 2 * (basis.order + 1); j++, jl++, j1++)
				A[il][jl] += Alocal.x[i1][j1];
		}

		for (size_t i = 0, il = 2 * x0 * basis.order; i < 2 * (basis.order + 1); i++, il++)
			b[il] += blocal[i];

		for (size_t i = 0, i1 = 2 * (basis.order + 1), il = 2 * x1 * basis.order; i < 2 * (basis.order + 1); i++, i1++, il++)
			b[il] += blocal[i1];
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
	void AddBorderInMatrix()
	{

		for (size_t i = 0; i < index.iyMax; i++)
		{
			int p = i * (2 * index.ixMax);
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[0].Us(x[0], y[i]);

			p++;
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[0].Uc(x[0], y[i]);
		}

		for (size_t i = 0; i < index.ixMax; i++)
		{
			int p = 2 * i;
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[1].Us(x[i], y[0]);

			p++;
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[1].Uc(x[i], y[0]);
		}

		for (size_t i = 0; i < index.iyMax; i++)
		{
			int p = 2 * ((i+1)* index.ixMax - 1);
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[2].Us(x[index.ixMax - 1], y[i]);

			p++;
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[2].Uc(x[index.ixMax - 1], y[i]);
		}


		for (size_t i = 0; i < index.ixMax; i++)
		{
			int p = 2 * (index.ixMax * (index.iyMax - 1) + i);
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[3].Us(x[i], y[index.iyMax - 1]);

			p++;
			for (size_t j = 0; j < m; j++)
				A[p][j] = 0;
			A[p][p] = 1;
			b[p] = border[3].Uc(x[i], y[index.iyMax - 1]);
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
