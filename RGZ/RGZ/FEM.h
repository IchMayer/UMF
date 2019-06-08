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
vector<double> operator/(vector<double> a, double b)
{
	for (size_t i = 0; i < a.size(); i++)
		a[i] /= b;
	return a;
}

//Разряженный формат (если нет, то плотный)
#define A1

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
		double b;
		function<double(double, double, double)>& U;

		//Задание первой или второй краевой функции
		Border(function<double(double, double, double)> &f, TypeBorder t) : U(f)
		{
			type = t;
			b = 1.0;
		}

		//Задание третьей краевой функции
		Border(function<double(double, double, double)> &f, TypeBorder t, double betta) : U(f)
		{
			type = t;
			b = betta;
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

		vector<double> localF(pow(basis.order + 1, 2));

		Matrix Glocal(9);
		Matrix Mlocal(9);

		double hx, hy;
		
		vector<double> ts(4);
		vector<double> th(4);
		
		vector<double> q1(9);	//значение q на k - 1 шаге по времени
		vector<double> q2(9);	//значение q на k - 2 шаге по времени
		vector<double> q3(9);	//значение q на k - 3 шаге по времени

		double mid_x, mid_y;
		Matrix Alocal;
		vector<double> localb;

		auto qlast = q;

		for (size_t k = u.size(); k < t.size(); k++)
		{
			ClearLastIter();
			lastIter = 0;
			q = u[u.size() - 1];
			qlast = q;

			ts = { ((t[k - 2] - t[k]) * (t[k - 1] - t[k])) / ((t[k - 3] - t[k - 2]) * (t[k - 3] - t[k - 1]) * (t[k - 3] - t[k])),
				((t[k - 3] - t[k]) * (t[k] - t[k - 1])) / ((t[k - 3] - t[k - 2])* (t[k - 2] - t[k - 1]) * (t[k - 2] - t[k])),
				((t[k - 3] - t[k]) * (t[k - 2] - t[k])) / ((t[k - 1] - t[k - 3]) * (t[k - 1] - t[k - 2]) * (t[k - 1] - t[k])),
				1.0 / (t[k] - t[k - 3]) + 1.0 / (t[k] - t[k - 2]) + 1.0 / (t[k] - t[k - 1])};

			th = { -2.0 * (t[k - 2] + t[k - 1] - 2.0 *t[k]) / ((t[k - 3] - t[k - 2]) * (t[k - 3] - t[k - 1]) * (t[k - 3] - t[k])),
				2.0 * (t[k - 3] + t[k - 1] - 2.0* t[k]) / ((t[k - 3] - t[k - 2])*(t[k - 2] - t[k - 1]) * (t[k - 2] - t[k])),
				2.0 * (t[k - 3] + t[k - 2] - 2.0*t[k]) / ((t[k - 3] - t[k - 1]) * (t[k - 1] - t[k - 2]) * (t[k - 1] - t[k])),
				2.0 * (t[k - 3] + t[k - 2] + t[k - 1] - 3.0 * t[k]) / ((t[k - 3] - t[k]) * (t[k] - t[k - 2]) * (t[k] - t[k - 1]))};

			do
			{
				for (size_t i = 0; i < x.size() - 1; i++)
				{
					for (size_t j = 0; j < y.size() - 1; j++)
					{
						hx = x[i + 1] - x[i];
						hy = y[j + 1] - y[j];

						for (size_t i = 0; i < 9; i++)
							for (size_t j = 0; j < 9; j++)
								Glocal.x[i][j] = (hy / hx * G.x[mu(i)][mu(j)] * M.x[nu(i)][nu(j)] + hx / hy * G.x[nu(i)][nu(j)] * M.x[mu(i)][mu(j)]);

						for (size_t i = 0; i < 9; i++)
							for (size_t j = 0; j < 9; j++)
								Mlocal.x[i][j] = hx * hy * M.x[mu(i)][mu(j)] * M.x[nu(i)][nu(j)];

						CreateLocalA(Glocal, Mlocal, Alocal, ts[3], th[3]);
						CreateLocalF(i, j, k, localF);
						CreateLocalQ(i, j, k, q1, q2, q3);

						CreateLocalb(Mlocal, localb, localF, q1, q2, q3, ts[2], ts[1], ts[0], th[2], th[1], th[0]);

						AddLocalMatrix(i, j, Alocal, localb);
					}
				}

				AddBorderInMatrix(k);

				//Prof2Raz();
				//high_resolution_clock::time_point t1 = high_resolution_clock::now();
				////lastIter = MSG(ig, jg, ggl, ggu, di, b, q1, eps);
				//lastIter = LOS(ig, jg, ggl, ggu, di, b, q1, eps);
				//high_resolution_clock::time_point t2 = high_resolution_clock::now();
				//tq1 = duration_cast<microseconds>(t2 - t1).count();

				//for (auto& i : A)
				//{
				//	for (auto& j : i)
				//		cout << j << " ";
				//	cout << endl;
				//}
				//cout << Norm(A * q - b) / Norm(b) << endl;


#ifdef A1
				MultAF(q, localb);
				if (Norm(localb - b) / Norm(b) < eps)
					break;
				qlast = q;
				//LOS(ig, jg, ggl, ggu, di, b, q, eps);
				//BSG_STAB(ig, jg, ggl, ggu, di, b, q, eps);
				MSG_LOS_BCGSTAB(ig, jg, ggl, ggu, di, b, q, eps);
#else
				if (Norm(A * q - b) / Norm(b) < eps)
					break;
				qlast = q;
				gause(A, b, q);
#endif // A1
			} while (Norm(q - qlast) / Norm(q) < eps);
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
	//Задать Сигму
	void SetSigma(double Sigma) { sigma = Sigma; }
	//Задать Лямбду
	void SetLambda(double l) { lambda = l; }
	//Задать Хи
	void SetHi(double Hi) { hi = Hi; }
	//Доьавление решения
	void AddSolution(vector<double> &q)
	{
		u.push_back(q);
	}
	//Добавление стороны
	void AddBorder(Border b)
	{
		border.push_back(b);
	}
	//Задание правой части
	void SetF(function<double(double, double, double)> func)
	{
		f = func;
	}
	//SetX
	void SetX(vector<double> X) { x = X; }
	//SetY
	void SetY(vector<double> Y) { y = Y; }
	//SetT
	void SetT(vector<double> T) { t = T; }
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

	vector<double>	x;				//Cетка по x
	vector<double>	y;				//Сетка по y
	vector<double>	t;				//Сетка по t

	double			lambda;			//Лямбда
	double			sigma;			//Сигма
	double			hi;				//Хи
	function<double(double, double, double)> f;

	//Результат
	vector<vector<double>>	u;		//Искомая функция

	//Переходные вычисления
	int m;

#ifndef A1
	vector<vector<double>> A;		//Матрица левой части в диагональном виде
#endif // !A1
	vector<double>	b;				//Вектора правой части				

	Matrix M;
	Matrix G;

	vector<double> di;	//Диагональные элементы матрицы А
	vector<double> ggu;	//Верхний треугольник матрицы А в разреженном формате
	vector<double> ggl;	//Нижний треугольний матрицы А в разреженном формате
	vector<size_t> ig;	// Массив индексов
	vector<size_t> jg;	// Другой массив индексов

	ShapeIndex		index;
	//Обработка введеных данных
	void Preprocessing()
	{
		if (t.size() < u.size())
			throw "Установлино неправильное число времен";

		if (u.size() < 3)
			throw "Введено недостаточно начальных условий";

		index.ixMax = 2 * x.size() - 1;
		index.iyMax = 2 * y.size() - 1;

		m = index.ixMax * index.iyMax;

		b.resize(m);
		q.resize(m);

		int order = basis.order + 1;

		di.resize(m);
		ig.resize(m + 1);
		
#ifdef A1
		for (size_t i = 0, i1 = 0; i < x.size() - 1; i++, i1 += 2)
			for (size_t j = 0, j1 = 0; j < y.size() - 1; j++, j1 += 2)
				for (size_t j2 = 0; j2 < 3; j2++)
					for (size_t i2 = 0; i2 < 3; i2++)
						for (size_t j3 = 0; j3 < 3; j3++)
							for (size_t i3 = 0; i3 < 3; i3++)
								CreateA1Elem((i1 + i3) + (j1 + j3) * index.ixMax, (i1 + i2) + (j1 + j2) * index.ixMax);
#else
		A.resize(m);
		for (size_t i = 0; i < m; i++)
			A[i].resize(m);
#endif // A1


		CreateMG();
	}
	//Создание матрицы массы и матрицы жесткости для 1D
	void CreateMG()
	{
			G.x = vector<vector<double>>{ {7, -8, 1}, {-8, 16, -8}, {1, -8, 7}};
			G = G  * (lambda /  3.0);
			G.n = 3;

			M.x = vector<vector<double>>{ {4, 2, -1}, {2, 16, 2}, {-1, 2, 4}};
			M = M / 30.0;
			M.n = 3;
	}
	//Создание локальной матрицы A
	inline void CreateLocalA(Matrix &G, Matrix &M, Matrix &Alocal, double ts, double th)
	{
		Alocal = G + M * (ts * sigma + th * hi);
	}
	//Создание локаольного вектора B
	inline void CreateLocalb(Matrix &M, vector<double>& blocal, vector<double> &localf,
						vector<double> &qlast1, vector<double> &qlast2, vector<double> &qlast3,
						double &ts1, double &ts2, double &ts3, double &th1, double &th2, double &th3)
	{
		blocal = M * (localf - (((qlast1 * ts1 + qlast2 * ts2 + qlast3 * ts3) * sigma) + ((qlast1 * th1 + qlast2 * th2 + qlast3 * th3) * hi)));
	}
	//Провкрка введенных данных
	void checkdata()
	{

	}
	//Добавление локальных матриц в глобальную
	void AddLocalMatrix(int i, int j, Matrix &localA, vector<double> &localb)
	{
		vector<int> local2full;
		for (size_t j1 = 0, j2 = 2 * j; j1 < 3; j1++, j2++)
			for (size_t i1 = 0, i2 = 2 * i; i1 < 3; i1++, i2++)
				local2full.push_back(i2 + j2 * index.ixMax);

		for (size_t j = 0; j < 9; j++)
		{
			for (size_t i = 0; i < 9; i++)
#ifdef A1
				AF(local2full[j], local2full[i]) += localA.x[j][i];
#else
				A[local2full[j]][local2full[i]] += localA.x[j][i];
#endif // A1
			b[local2full[j]] += localb[j];
		}
	}
	//Очистка A и b c предыдущей итерации
	void ClearLastIter()
	{
#ifdef A1
		for (size_t i = 0; i < ggu.size(); i++)
		{
			ggu[i] = 0;
			ggl[i] = 0;
		}
		for (auto& i : di)
			i = 0;
#else
		for(auto& i : A)
			for (auto& j : i)
				j = 0;
#endif // A1
		for (size_t i = 0; i < b.size(); i++)
			b[i] = 0;
	}
	//Очищение места где должны быть краевые условия
	void SetThirdBorderA(double hx, int p1, int p2, int p3)
	{
#ifdef A1
		AF(p1,p1) += M.x[0][0] * hx; AF(p1,p2) += M.x[0][1] * hx; AF(p1,p3) += M.x[0][2] * hx;
		AF(p2,p1) += M.x[1][0] * hx; AF(p2,p2) += M.x[1][1] * hx; AF(p2,p3) += M.x[1][2] * hx;
		AF(p3,p1) += M.x[2][0] * hx; AF(p3,p2) += M.x[2][1] * hx; AF(p3,p3) += M.x[2][2] * hx;
#else
		A[p1][p1] += M.x[0][0] * hx; A[p1][p2] += M.x[0][1] * hx; A[p1][p3] += M.x[0][2] * hx;
		A[p2][p1] += M.x[1][0] * hx; A[p2][p2] += M.x[1][1] * hx; A[p2][p3] += M.x[1][2] * hx;
		A[p3][p1] += M.x[2][0] * hx; A[p3][p2] += M.x[2][1] * hx; A[p3][p3] += M.x[2][2] * hx;
#endif // A1
	}
	//Добавление краевых условаий в матрицу А и b
	void AddBorderInMatrix(int k)
	{
		//ClearBorder();

		int p;
		int p1;
		int p2;
		//Левая сторона не 1 краевое
		if (border[0].type != border[0].First)
			for (size_t i = 0; i < y.size() - 1; i++)
			{
				p = (2 * i) * index.ixMax;
				p1 = (2 * i + 1) * index.ixMax;
				p2 = (2 * i + 2) * index.ixMax;

				double hy = y[i + 1] - y[i];
				if (border[0].type == border[0].Third)
				{
					hy *= border[0].b;
					SetThirdBorderA(hy, p, p1, p2);
				}

				vector<double> loc2b = M * vector<double>{	border[0].U(x[0], y[i], t[k]),
															border[0].U(x[0], (y[i + 1] + y[i]) / 2.0, t[k]),
															border[0].U(x[0], y[i + 1], t[k])};
				b[p] += loc2b[0] * hy;
				b[p1] += loc2b[1] * hy;
				b[p2] += loc2b[2] * hy;

			}

		//Нижняя сторона не 1 краевое
		if (border[1].type != border[1].First)
			for (size_t i = 0; i < x.size() - 1; i++)
			{
				p = 2 * i;
				p1 = 2 * i + 1;
				p2 = 2 * i + 2;

				double hx = x[i + 1] - x[i];
				if (border[1].type == border[1].Third)
				{
					hx *= border[1].b;
					SetThirdBorderA(hx, p, p1, p2);
				}

				vector<double> loc2b = M * vector<double>{	border[1].U(x[i], y[0], t[k]),
															border[1].U((x[i + 1] + x[i]) / 2, y[0], t[k]),
															border[1].U(x[i + 1], y[0], t[k])};

				b[p] += loc2b[0]  * hx;
				b[p1] += loc2b[1] * hx;
				b[p2] += loc2b[2] * hx;

			}

		//Правая сторона не 1 краевое
		if (border[2].type != border[2].First)
			for (size_t i = 0; i < y.size() - 1; i++)
			{
				p = ((2 * i + 1)* index.ixMax - 1);
				p1 = ((2 * i + 2)* index.ixMax - 1);
				p2 = ((2 * i + 3)* index.ixMax - 1);

				double hy = y[i + 1] - y[i];
				if (border[2].type == border[2].Third)
				{
					hy *= border[2].b;
					SetThirdBorderA(hy, p, p1, p2);
				}

				vector<double> loc2b = M * vector<double>{	border[2].U(x[x.size() - 1], y[i], t[k]),
															border[2].U(x[x.size() - 1], (y[i + 1] + y[i]) / 2.0, t[k]),
															border[2].U(x[x.size() - 1], y[i + 1], t[k])};

				b[p] += loc2b[0]  * hy;
				b[p1] += loc2b[1] * hy;
				b[p2] += loc2b[2] * hy;

			}

		//Верхняя сторона не 1 краевое
		if (border[3].type != border[3].First)
			for (size_t i = 0; i < x.size() - 1; i++)
			{
				p = index.ixMax * (index.iyMax - 1) + 2 * i;
				p1 = index.ixMax * (index.iyMax - 1) + 2 * i + 1;
				p2 = index.ixMax * (index.iyMax - 1) + 2 * i + 2;
				
				double hx = x[i + 1] - x[i];
				if (border[3].type == border[3].Third)
				{
					hx *= border[3].b;
					SetThirdBorderA(hx, p, p1, p2);
				}

				vector<double> loc2b = M * vector<double>{	border[3].U(x[i], y[y.size() - 1], t[k]),
															border[3].U((x[i + 1] + x[i]) / 2, y[y.size() - 1], t[k]),
															border[3].U(x[i + 1], y[y.size() - 1], t[k])};

				b[p] += loc2b[0] * hx;
				b[p1] += loc2b[1] * hx;
				b[p2] += loc2b[2] * hx;

			}

		//Если левая сторона 1 краевое
		if(border[0].type == border[0].First)
			for (size_t i = 0; i < index.iyMax; i++)
			{
				p = i * index.ixMax;
				for (size_t j = 0; j < m; j++)
#ifdef A1
					AF(p, j, 0);
				AF(p, p) = 1;
#else
					A[p][j] = 0;
				A[p][p] = 1;
#endif // A1

				if (i % 2 == 0)
					b[p] = border[0].U(x[0], y[i / 2], t[k]);
				else
					b[p] = border[0].U(x[0], (y[i / 2] + y[i / 2 + 1]) / 2, t[k]);
			}

		//Если нижняя сторона 1 краевое
		if (border[1].type == border[1].First)
			for (size_t i = 0; i < index.ixMax; i++)
			{
				p = i;
				for (size_t j = 0; j < m; j++)
#ifdef A1
					AF(p, j, 0);
				AF(p, p) = 1;
#else
					A[p][j] = 0;
				A[p][p] = 1;
#endif // A1

				if (i % 2 == 0)
					b[p] = border[1].U(x[i / 2], y[0], t[k]);
				else
					b[p] = border[1].U((x[i / 2] + x[i / 2 + 1]) / 2, y[0], t[k]);
			}

		//Если правая сторона 1 краевое
		if (border[2].type == border[2].First)
			for (size_t i = 0; i < index.iyMax; i++)
			{
				p = ((i + 1)* index.ixMax - 1);
				for (size_t j = 0; j < m; j++)
#ifdef A1
					AF(p, j, 0);
				AF(p, p) = 1;
#else
					A[p][j] = 0;
				A[p][p] = 1;
#endif // A1

				if (i % 2 == 0)
					b[p] = border[2].U(x[(index.ixMax - 1) / 2], y[i / 2], t[k]);
				else
					b[p] = border[2].U(x[(index.ixMax - 1) / 2], (y[i / 2] + y[i / 2 + 1]) / 2, t[k]);
			}

		//Если верхняя сторона 1 краевое
		if (border[3].type == border[3].First)
			for (size_t i = 0; i < index.ixMax; i++)
			{
				p = index.ixMax * (index.iyMax - 1) + i;
				for (size_t j = 0; j < m; j++)
#ifdef A1
					AF(p, j, 0);
				AF(p, p) = 1;
#else
					A[p][j] = 0;
				A[p][p] = 1;
#endif // A1
				if (i % 2 == 0)
					b[p] = border[3].U(x[i / 2], y[(index.iyMax - 1) / 2], t[k]);
				else
					b[p] = border[3].U((x[i / 2] + x[i / 2 + 1]) / 2, y[(index.iyMax - 1) / 2], t[k]);
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

	//Перенос в 2D матрицу G и M из 1D
	inline int mu(int i){return (i % 3);}
	inline int nu(int i) {	return i / 3;}

	//Создаем на каждом шаге локальный вектор правой части
	void CreateLocalF(int i, int j, int k, std::vector<double>& localF)
	{
		double mid_x = (x[i + 1] + x[i]) / 2.0;
		double mid_y = (y[j + 1] + y[j]) / 2.0;

		localF[0] = f(x[i],		y[j], t[k]);
		localF[1] = f(mid_x,	y[j], t[k]);
		localF[2] = f(x[i + 1], y[j], t[k]);

		localF[3] = f(x[i],		mid_y, t[k]);
		localF[4] = f(mid_x,	mid_y, t[k]);
		localF[5] = f(x[i + 1], mid_y, t[k]);

		localF[6] = f(x[i],		y[j + 1], t[k]);
		localF[7] = f(mid_x,	y[j + 1], t[k]);
		localF[8] = f(x[i + 1], y[j + 1], t[k]);
	}
	//Создание на каждом шаге локальный вектор решений на предыдущих итерациях
	void CreateLocalQ(int i, int j, int k, std::vector<double>& q1, std::vector<double>& q2, std::vector<double>& q3)
	{
		for (size_t j1 = 0, j2 = 2 * j, p = 0; j1 < 3; j1++, j2++)
			for (size_t i1 = 0, i2 = 2 * i; i1 < 3; i1++, i2++, p++)
			{
				q1[p] = u[k - 1][i2 + j2 * index.ixMax];
				q2[p] = u[k - 2][i2 + j2 * index.ixMax];
				q3[p] = u[k - 3][i2 + j2 * index.ixMax];
			}
	}


#ifdef A1
	//Добавлеие элемента в разряженную матрицу 
	void CreateA1Elem(int i, int j)
	{
		if (i == j)
			return;

		if (j > i)
		{
			int c = j;
			j = i;
			i = c;
		}

		for (size_t k = ig[i]; k < ig[i + 1]; k++)
			if (jg[k] == j)
				return;

		for (size_t i1 = i + 1; i1 < ig.size(); i1++)
			ig[i1]++;

		jg.push_back(0);
		for (int k = ig[ig.size() - 1] - 1; k >= ig[i + 1]; k--)
			jg[k] = jg[k - 1];

		ggu.push_back(0);
		ggl.push_back(0);

		for (int k = ig[i]; k < ig[i + 1] - 1; k++)
		{
			if (jg[k] > j)
			{
				for (size_t p = ig[i + 1] - 1; p > k; p--)
					jg[p] = jg[p - 1];
				jg[k] = j;
				return;
			}
		}
		jg[ig[i + 1] - 1] = j;
	}

	//Ображение к элементку в разряженной матрице
	double& AF(int i, int j)
	{
		if (i == j)
			return di[i];
		if (i > j)
		{
			for (int k = ig[i]; k < ig[i + 1]; k++)
			{
				if (jg[k] == j)
					return ggl[k];
			}
		}
		else
		{
			for (int k = ig[j]; k < ig[j + 1]; k++)
			{
				if (jg[k] == i)
					return ggu[k];
			}
		}
		double p = 0;
		return p;
	}

	//Умножение матрицы в разряженном формате на вектор
	void MultAF(vector<double> &b, vector<double> &result)
	{
		result.resize(m);
		
		for (int i = 0; i < m; i++)
		{
			result[i] = di[i] * b[i];
			for (int k = ig[i]; k < ig[i + 1]; k++)
			{
				int j = jg[k];
				result[i] += ggl[k] * b[j];
				result[j] += ggu[k] * b[i];
			}
		}
	}

	//Задать значение элементу матрицы
	void AF(int i, int j, double c)
	{
		if (i == j)
		{
			di[i] = c;
			return;
		}
		if (i > j)
		{
			for (int k = ig[i]; k < ig[i + 1]; k++)
			{
				if (jg[k] == j)
				{
					ggl[k] = c;
					return;
				}
			}
		}
		else
		{
			for (int k = ig[j]; k < ig[j + 1]; k++)
			{
				if (jg[k] == i)
				{
					ggu[k] = c;
					return;
				}
			}
		}
	}
#endif // A1
};
