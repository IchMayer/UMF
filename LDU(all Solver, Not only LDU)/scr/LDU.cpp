// LDU.cpp : Определяет экспортированные функции для приложения DLL.
//

#include "stdafx.h"


namespace Matrix
{
	enum class Errors
	{
		noError,
		ErrorOpenFile,
		MatrixNotCount
	};

	enum class Format
	{
		LDU,
		Profile
	};

	template<class T>
	class Matrix
	{
	public:
		Matrix() {}
		~Matrix() {}

		std::vector<T> al, au, di, f, x;
		std::vector<size_t> ia;

		size_t n, m;
		Format format;
		Errors error;
	private:

	};

	namespace Gause
	{
		template<class T>
		class Gause : Matrix<T>
		{
		public: 
			std::vector<std::vector<T>> A;
			std::vector<double> f, x;

			Gause()
			{
				this->format = Format::Profile;
				error = Errors::noError;
			}
			~Gause() {}
			Errors error;

		private:
			size_t n; //Размерность матрицы
		};

		template<class T>
		void GausBack(Gause<T> &matrix)
		{
			std::vector<T> &X = matrix.x;
			std::vector<T> &F = matrix.f;
			std::vector<std::vector<T>> &A = matrix.A;
			auto n = F.size();
			X.resize(n);
			if (matrix.error != Errors::noError)
				return;
			for (int i = n - 1; i >= 0; i--)
			{
				double sum = 0;
				for (size_t j = i + 1; j < n; j++)
					sum += A[i][j] * X[j];
				X[i] = (F[i] - sum) / A[i][i];
			}
		}

		template<class T>
		void CountUpTrianMatrix(Gause<T> &matrix)
		{
			if (matrix.error != Errors::noError)
				return;
			std::vector<std::vector<T>> &A = matrix.A;
			std::vector<T> &F = matrix.f;
			auto n = F.size();
			for (size_t i = 0; i < n; i++)
			{
				//Опредиление главного Элемента в столбце
				T Max = A[i][i];
				size_t index = i;
				for (size_t j = i + 1; j < n; j++)
				{
					if (abs(A[j][i]) > abs(Max))
					{
						Max = A[j][i];
						index = j;
					}
				}

				//Проверка единственность
				if (!Max)
				{
					matrix.error = Errors::MatrixNotCount;
					return;
				}

				//Свап строки
				if (i != index)
				{
					A[index].swap(A[i]);
					T buf = F[i];
					F[i] = F[index];
					F[index] = buf;
				}

				//Постройка верхней треугольной матрицы
				for (size_t j = i + 1; j < n; j++)
				{
					if (!A[j][i])
						continue;
					T c = A[j][i] / Max;
					A[j][i] = 0.;
					for (size_t l = i + 1; l < n; l++)
						A[j][l] -= A[i][l] * c;
					F[j] -= F[i] * c;
				}
			}
		}
		
		template<class T>
		Gause<T> Create(std::vector<double> &al, std::vector<double> &au, std::vector<double> &di, std::vector<size_t> &ia)
		{
			Gause<T> m;
			std::vector<std::vector<T>> &A = m.A;
			if (!di.size())
			{
			 	m.Error = Errors::MatrixNotCount;
				return;
			}
			else
			{
				A.resize(m.n);
				for (size_t i = 0; i < m.n; A[i++].resize(m.n));
			}
			
			m.n = di.size();
			auto n = m.n;
			for (size_t i = 0; i < n; i++)
				A[i][i] = di[i];

			for (size_t i = 0; i < n; i++)
			{
				for (size_t m = ia[i], j = i - ia[i + 1] + m; m < ia[i + 1]; m++, j++)
				{
					A[i][j] = al[m];
					A[j][i] = au[m];
				}
			}
		}

	}

	namespace LDU
	{
		template<class T>
		class LDU: Matrix<T>
		{
		public:
			LDU<T>()
			{
				this->format = Format::LDU;
			}
			~LDU<T>() {}
			std::vector<T> al, au, di, f, x;
			std::vector<size_t> ia;

			size_t n, m;
			Format format;
			Errors error;
		};


		template<class T>
		void CountX(LDU<T> &matrix, std::vector<T> F)
		{
			size_t n = matrix.n;
			std::vector<T> &D = matrix.di;
			std::vector<T> &U = matrix.au;
			std::vector<T> &L = matrix.al;
			std::vector<T> &X = matrix.x;
			Errors Error = matrix.error;
			std::vector<size_t> &ia = matrix.ia;
			for (auto &elem : X)
				elem = 0;

			if (Error != Errors::noError)
				return;
			//Поиск Y
			std::vector<T> &Y = X;
			for (size_t i = 0; i < n; i++)
			{
				double sum = 0;
				for (size_t m = ia[i],
					j = i - ia[i + 1] + m;
					m < ia[i + 1]; m++)
					sum += Y[j++] * L[m];
				Y[i] += F[i] - sum;
			}

			//Поиск Z
			std::vector<T> &Z = X;
			for (size_t i = 0; i < n; i++)
				Z[i] /= D[i];

			//Поиск X
			for (int i = n - 1; i >= 0; i--)
				for (int j = ia[i + 1] - ia[i] - 1,
					j1 = i - j - 1,
					j2 = ia[i + 1] - j - 1;
					j >= 0; j--)
					X[j1++] -= X[i] * U[j2++];
		}

		template<class T>
		void CountLDU(LDU<T> &matrix)
		{
			size_t n = matrix.n;
			size_t m = matrix.m;
			auto &D = matrix.di;
			auto &U = matrix.au;
			auto &L = matrix.al;
			auto &Error = matrix.error;
			auto &ia = matrix.ia;

			double CompareNumber;
			if (std::is_same<T, double>::value)
				CompareNumber = pow(10, 15);
			else
				CompareNumber = pow(10, 7);

			for (size_t i = 0; i < n; i++)
			{
				double sum = 0;
				for (size_t j = ia[i], j1 = i - ia[i + 1] + j; j < ia[i + 1]; j++)
				{
					double sum1 = 0, sum2 = 0;
					for (size_t k = min(j - ia[i], ia[i - ia[i + 1] + j + 1] - ia[i - ia[i + 1] + j]),
						i1 = j - k,
						i2 = i - ia[i + 1] + j - k,
						i3 = ia[i - ia[i + 1] + j + 1] - k;
						k > 0; k--)
					{
						sum1 += L[i1] * D[i2] * U[i3];
						sum2 += L[i3++] * D[i2++] * U[i1++];
					}
					L[j] = (L[j] - sum1) / D[j1];
					U[j] = (U[j] - sum2) / D[j1];
					sum += L[j] * D[j1++] * U[j];
				}
				//Если сделать проверка на inf после деления на 0, то можем получить не то что нам нужно
				//Т.к. при разность (D[i] - sum) может получиться не 0, а число близкое к нему(разность нецелого типа)
				//и в таком случае при делении числа на эту разность может получиться не inf.
				if (abs((D[i] - sum)) < abs(D[i] / CompareNumber))
				{
					Error = Errors::MatrixNotCount;
					return;
				}
				D[i] -= sum;
			}
		}

		template <class T>
		void CountLDU(LDU<T> &matrix, std::vector<T> &al, std::vector<T> &au, std::vector<T> &di, std::vector<size_t> &ia)
		{
			size_t n = matrix.n;
			size_t m = matrix.m;
			auto &D = matrix.di;
			auto &U = matrix.au;
			auto &L = matrix.al;
			D.resize(n);
			U.resize(m);
			L.resize(m);
			auto &Error = matrix.error;
			
			matrix.ia = ia;

			double CompareNumber;
			if (std::is_same<T, double>::value)
				CompareNumber = pow(10, 15);
			else
				CompareNumber = pow(10, 7);

			for (size_t i = 0; i < n; i++)
			{
				double sum = 0;
				for (size_t j = ia[i], j1 = i - ia[i + 1] + j; j < ia[i + 1]; j++)
				{
					double sum1 = 0, sum2 = 0;
					for (size_t k = min(j - ia[i], ia[i - ia[i + 1] + j + 1] - ia[i - ia[i + 1] + j]),
						i1 = j - k,
						i2 = i - ia[i + 1] + j - k,
						i3 = ia[i - ia[i + 1] + j + 1] - k;
						k > 0; k--)
					{
						sum1 += L[i1] * D[i2] * U[i3];
						sum2 += L[i3++] * D[i2++] * U[i1++];
					}
					L[j] = (al[j] - sum1) / D[j1];
					U[j] = (au[j] - sum2) / D[j1];
					sum += L[j] * D[j1++] * U[j];
				}
				//Если сделать проверка на inf после деления на 0, то можем получить не то что нам нужно
				//Т.к. при разность (D[i] - sum) может получиться не 0, а число близкое к нему(разность нецелого типа)
				//и в таком случае при делении числа на эту разность может получиться не inf.
				if (abs((D[i] - sum)) < abs(D[i] / CompareNumber))
				{
					Error = Errors::MatrixNotCount;
					return;
				}
				D[i] = di[i] - sum;
			}
		}

		template<class U>
		bool Open(std::string path, std::vector<U> &a, size_t len)
		{
			std::ifstream ia;
			ia.exceptions(std::ifstream::failbit | std::ifstream::badbit);
			try
			{
				ia.open(path);
				a.resize(len);
				for (size_t i = 0; i < len; i++)
					ia >> a[i];
				ia.close();
			}
			catch (std::ifstream::failure e)
			{
				std::cout << "Error Open File: " << path << std::endl;
				return false;
			}
			return true;
		}

		template<class T>
		void OpenFile(LDU<T> &matrix, std::string path)
		{
			int error = 0;
			size_t n, m;
			if (path.size())
				path += "/";
			std::vector<size_t> info(2);
			error += 
				!Open(path + "info.txt", info, 2);
			n = info[0], m = info[1];
			matrix.x.resize(n);
			error +=
				!Open(path + "al.txt", matrix.al, m) +
				!Open(path + "au.txt", matrix.au, m) +
				!Open(path + "di.txt", matrix.di, n) +
				!Open(path + "ia.txt", matrix.ia, n + 1) +
				!Open(path + "F.txt", matrix.f, n);
			if(error)
				matrix.error = Errors::ErrorOpenFile;
		}

		template<class T>
		LDU<T> Create(std::string path)
		{
			LDU<T> m;
			OpenFile(m, path);
			CountLDU(m);
			return m;
		}

		template<class T>
		LDU<T> Create(std::vector<double> &al, std::vector<double> &au, std::vector<double> &di, std::vector<size_t> &ia)
		{
			LDU<T> m;
			CountLDU(m, al, au, di, ia);
			return m;
		}
	}

	namespace MSG_LOS_BCGSTAB
	{
		template<class T>
		class Matrix
		{
		public:
			Matrix() { Maxiter = 10000; }
			~Matrix() {}

			void Init(size_t n, size_t m, std::vector<size_t> iig, std::vector<size_t> ijg, std::vector<double> gl, std::vector<double> gu, std::vector<double> Di, std::vector<double> F, double EPS)
			{
				N = n;
				M = m;

				e = EPS;

				ig = iig;
				jg = ijg;

				result.resize(N);

				ggl = gl;
				ggu = gu;
				
				di = Di;

				f = F;
			}

			//Set matrix dimension
			///<param name = "n"> New matrix dimension</param>
			void SetN(size_t n) { N = n; }
			//Get matrix dimention
			size_t n() { return N; }

			//Возращает резултат операций
			std::vector<T> Result() { return result; }
			//Возращает количество операций
			size_t K() { return k; }

			//Set max iteration
			///<param name = "Max"> New max ieration </param>
			void Setmaxiter(size_t Max) { Maxiter = Max; }
			//Get max iteration
			size_t maxiter() { return Maxiter; }

			///<param name = "i"> i = 0 без предобуславливония
			///i = 1 диагональное
			///i = 2 Холесского</param>
			void preconditioning(int i)
			{
				switch (i)
				{
					//Без предобуславливония
				case 0:

					break;
					//Диагональное
				case 1:
					//			factorizationLU(false);
					factorizationLLT();
					break;
					//Холесского
			//	case 2:
				//	factorizationLLT(true);
					//break;
					//LU
				case 3:
					factorizationLU(true);
					break;
				default:
					break;
				}
			}

			//i = 0 без предобуславливония
			//i = 1 диагональное
			//i = 2 Холесского</param>
			void preconditioning()
			{
				preconditioning(3);
			}

			void BSG_STAB(int metod)
			{
				NormF = Norm(f);
				switch (metod)
				{
				case 0:
				{
					std::vector<T> Best = result;
					r = Residual(f, MultMatrixVector(Best));
					auto p = r;
					z = r;
					auto s = r;
					std::vector<T> d;

					T rpScolar = Scolar(r, p);
					T rpLastScolar;

					T alfa, betta;
					T CoudA;
					T min = 0.;
					while (checkEnd(rpScolar, CoudA))
					{
						d = MultMatrixVector(z);

						alfa = rpScolar / Scolar(s, d);

						Best = Sum(Best, MultVectorOnT(z, alfa));

						r = Residual(r, MultVectorOnT(d, alfa));

						d = MultTMatrixVector(s);
						p = Residual(p, MultVectorOnT(d, alfa));

						rpLastScolar = rpScolar;
						rpScolar = Scolar(r, p);

						betta = rpScolar / rpLastScolar;

						z = Sum(r, MultVectorOnT(z, betta));
						z = Sum(p, MultVectorOnT(s, betta));

						k++;
						if (min > CoudA || !min)
						{
							min = CoudA;
							result = Best;
						}
						k++;
					}

					break;
				}
				case 1:
				case 2:
					break;
				case 3:
				{
					auto Best = result;
					x = Up(U, D, result);
					r = Residual(f, MultMatrixVector(x));
					r = Down(L, r);
					z = r;
					auto p = r;
					auto s = r;

					T rpScolar = Scolar(r, p);
					T rpLastScolar;

					T alfa;
					T betta;
					NormF = Norm(f);

					std::vector<T> d;
					std::vector<T> buf;
					std::vector<T> buf1;
					T Coud;
					T min = 0.;
					while (checkEnd(rpScolar, Coud))
					{
						d = Up(U, D, z);
						d = MultMatrixVector(d);
						d = Down(L, d);

						double c = Scolar(s, d);
						
						if (!c)
							break;

						alfa = rpScolar / c;

						buf1 = MultVectorOnT(z, alfa);
						x = Sum(x, buf1);

						buf1 = MultVectorOnT(d, alfa);
						r = Residual(r, buf1);

						buf = Up(L, s);
						buf = MultTMatrixVector(buf);
						buf = Down(U, D, buf);
						
						buf1 = MultVectorOnT(buf, alfa);
						p = Residual(p, buf1);
						
						rpLastScolar = rpScolar;
						rpScolar = Scolar(p, r);

						betta = rpScolar / rpLastScolar;

						buf1 = MultVectorOnT(z, betta);
						z = Sum(r, buf1);
						buf1 = MultVectorOnT(s, betta);
						s = Sum(p, buf1);

						if (min > Coud || !min)
						{
							min = Coud;
							Best = x;
						}

						k++;
					}

					result = Up(U, D, Best);

					//for (auto &i : result)
					//	i = 0;

					//for (size_t i = 0; i < N; i++)
					//{
					//	for (size_t j = ig[i]; j < ig[i + 1]; j++)
					//		result[jg[j]] += U[j] * Best[i];
					//	result[i] += D[i] * Best[i];
					//}
					break;
				}
				default:
					break;
				}
			}

			void ConjugateGradient(int metod)
			{
				NormF = Norm(f);
				switch (metod)
				{
					//Без
				case 0:
				{
					std::vector<T> Best = result;
					f = MultTMatrixVector(f);
					r = Residual(f, MultTMatrixVector(MultMatrixVector(Best)));
					z = r;
					T rLastScolar = Scolar(r);
					T CoudA;
					T min = 0.;
					while (checkEnd(rLastScolar, CoudA))
					{
						std::vector<T> d = MultTMatrixVector(MultMatrixVector(z));
						a = rLastScolar / Scolar(d, z);
						Best = Sum(Best, MultVectorOnT(z, a));
						rLast = r;
						r = Residual(r, MultVectorOnT(d, a));
						T lol = Scolar(r);
						b = lol / rLastScolar;
						rLastScolar = lol;
						z = Sum(r, MultVectorOnT(z, b));
						if (min > CoudA || !min)
						{
							min = CoudA;
							result = Best;
						}
						k++;
					}
					break;
				}
				// Диаганальное
				case 1:
					//LLT
				case 2:
				{
					std::vector<T> d = Residual(f, MultMatrixVector(result)); //r = U^-t * A ^t * L^ -T * L-1(f - Ax)
					d = Down(L, D, d);// L d1 = d
					d = Up(L, D, d);	//L^t d1 = d
					d = MultTMatrixVector(d); // A^t * d = d1
					r = Down(L, D, d);	// U^t r = d
					z = r;

					x.resize(N);
#pragma region Умножение_x-_=_Ux
					for (size_t i = 0; i < N; i++)
					{
						for (size_t j = ig[i]; j < ig[i + 1]; j++)
							x[i] += L[j] * result[jg[j]];
						x[i] += D[i] * result[i];
					}
#pragma endregion

					k = 0;
					//xbest = x;

					std::vector<T> Best = x;
					T rLastScolar = Scolar(r);
					T CoudA;
					T min = 0.;

					while (checkEnd(rLastScolar, CoudA))
					{

						std::vector<T> d = Up(L, D, z);
						d = MultMatrixVector(d);
						d = Down(L, D, d);
						d = Up(L, D, d);
						d = MultTMatrixVector(d);
						d = Down(L, D, d);
						T s = Scolar(d, z);
						if (!s)
							break;
						a = rLastScolar / s;
						x = Sum(x, MultVectorOnT(z, a));
						rLast = r;
						r = Residual(r, MultVectorOnT(d, a));
						T rScolar = Scolar(r);
						b = rScolar / rLastScolar;
						z = Sum(r, MultVectorOnT(z, b));
						rLastScolar = rScolar;
						k++;
						if (min > CoudA || !min)
						{
							min = CoudA;
							Best = x;
						}
					}
					result = Up(L, D, Best);
					break;
				}
				//LU
				case 3:
				{
					std::vector<T> d = Residual(f, MultMatrixVector(result)); //r = U^-t * A ^t * L^ -T * L-1(f - Ax)
					d = Down(L, d);// L d1 = d
					d = Up(L, d);	//L^t d1 = d
					d = MultTMatrixVector(d); // A^t * d = d1
					r = Down(U, D, d);	// U^t r = d
					z = r;

					x.resize(N);
#pragma region Умножение_x-_=_Ux
					for (size_t i = 0; i < N; i++)
					{
						for (size_t j = ig[i]; j < ig[i + 1]; j++)
							x[i] += U[j] * result[jg[j]];
						x[i] += D[i] * result[i];
					}
#pragma endregion

					k = 0;
					xbest = x;

					std::vector<T> Best = x;
					T rLastScolar = Scolar(r);
					T Coud;
					T min = 0.;
					while (checkEnd(rLastScolar, Coud))
					{

						std::vector<T> d = Up(U, D, z);
						d = MultMatrixVector(d);
						d = Down(L, d);
						d = Up(L, d);
						d = MultTMatrixVector(d);
						d = Down(U, D, d);
						T s = Scolar(d, z);
						if (!s)
							break;
						a = rLastScolar / s;
						x = Sum(x, MultVectorOnT(z, a));
						r = Residual(r, MultVectorOnT(d, a));
						T rScolar = Scolar(r);
						b = rScolar / rLastScolar;
						rLastScolar = rScolar;
						z = Sum(r, MultVectorOnT(z, b));

						if (min > Coud || !min)
						{
							min = Coud;
							Best = x;
						}

						k++;
					}
					result = Up(U, D, Best);
					break;
				}
				default:
					break;
				}
			}
			void LocalOptimalScheme(int metod)
			{
				NormF = Norm(f);
				switch (metod)
				{
				case 0:
				{
					r = Residual(f, MultMatrixVector(result));
					z = r;
					p = MultMatrixVector(z);
					k = 0;
					T ScolarP = Scolar(p);
					T ScolarR = Scolar(r);
					//while (checkEnd(ScolarR, ScolarP, ScolarR))
					while (checkEnd(r))
					{
						if (!ScolarP)
							break;
						a = Scolar(p, r) / ScolarP;
						result = Sum(result, MultVectorOnT(z, a));
						r = Residual(r, MultVectorOnT(p, a));
						T gg = Scolar(r);
						std::vector<T> d = MultMatrixVector(r);
						b = -1.0 * Scolar(p, d) / ScolarP;
						z = Sum(r, MultVectorOnT(z, b));
						p = Sum(d, MultVectorOnT(p, b));
						ScolarP = Scolar(p);
						k++;
					}
					break;
				}
				//DLLT
				case 1:
					//LLT
				case 2:
				{
					r = Down(L, D, Residual(f, MultMatrixVector(result)));
					z = Up(L, D, r);
					p = Down(L, D, MultMatrixVector(z));

					k = 0;
					while (checkEnd(r))
					{
						T s = Scolar(p);
						if (!s)
							break;
						a = Scolar(p, r) / s;
						result = Sum(result, MultVectorOnT(z, a));
						r = Residual(r, MultVectorOnT(p, a));
						std::vector<T> d = Up(L, D, r);
						d = MultMatrixVector(d);
						d = Down(L, D, d);
						b = -1.0 * Scolar(p, d) / s;
						z = Sum(Up(L, D, r), MultVectorOnT(z, b));
						p = Sum(d, MultVectorOnT(p, b));
						k++;
					}
					break;
				}
				//LU
				case 3:
				{
					r = Down(L, Residual(f, MultMatrixVector(result)));
					z = Up(U, D, r);
					p = Down(L, MultMatrixVector(z));

					k = 0;
					while (checkEnd(r))
					{
						T s = Scolar(p);
						if (!s)
							break;
						a = Scolar(p, r) / s;
						result = Sum(result, MultVectorOnT(z, a));
						r = Residual(r, MultVectorOnT(p, a));
						std::vector<T> d = Up(U, D, r);
						d = MultMatrixVector(d);
						d = Down(L, d);
						b = -1.0 * Scolar(p, d) / s;
						z = Sum(Up(U, D, r), MultVectorOnT(z, b));
						p = Sum(d, MultVectorOnT(p, b));
						k++;
					}
					break;
				}
				default:
					break;
				}
			}

		private:
			std::vector<T> f; // Массив правой части
			std::vector<T> di;	//Диагональные элементы матрицы А
			std::vector<T> ggu;	//Верхний треугольник матрицы А в разреженном формате
			std::vector<T> ggl;	//Нижний треугольний матрицы А в разреженном формате
			std::vector<T> L; // Матрица L
			std::vector<T> U; // НЕ диагональыне элементы  матрицы U
			std::vector<T> D; // Диагональные элементы матрицы U
			std::vector<size_t> ig;// Массив индексов
			std::vector<size_t> jg; // Другой массив индексов

			std::vector<T> result;	//Результат x
			std::vector<T> r, rLast, z, p, x, xbest;

			size_t N, M, Maxiter, k;
			T e, a, b, NormF;

			//Проверка окончания по невязки и по maxiter

			///<param name = "ScolarR"> Входной параметр. Сколяр вектора R</param>
			///<param name = "Return"> Выходной параметр. Невязка</param>
			//Проверка конца у МСГ
			bool checkEnd(T ScolarR, T &Return)
			{
				if (k == Maxiter)
					return false;
				Return = sqrt(ScolarR) / NormF;
				if (e > Return)
					return false;
				return true;
			}


			///<param name = "ScolarR"> Входной параметр. Сколяр вектора R</param>
			///<param name = "ScolarP"> Входной параметр. Сколяр вектора P</param>
			///<param name = "ReturnR"> Выходной параметр. Сколяр вектора R</param>
			//Проверка конца у ЛОС
			bool checkEnd(T ScolarR, T ScolarP, T &ReturnR)
			{
				if (k == Maxiter)
					return false;

				T R = sqrt(ScolarR) / NormF;
				T g = Scolar(r);
				if (e > R)
					return false;
				a = Scolar(p, r) / Scolar(p);
				T h = pow(a, 2) * ScolarP;
				ReturnR = ScolarR - h;
				return true;
			}

			bool checkEnd(std::vector<T> r)
			{
				if (k == Maxiter)
					return false;
				if (e > Norm(r) / NormF)
					return false;
				return true;
			}

			//Прямой ход
			std::vector<T> Down(std::vector<T> Matrix, std::vector<T> Diagonal, std::vector<T> R)
			{
				std::vector<T> x(N);
				for (size_t i = 0; i < N; i++)
				{
					T sum = 0;
					for (size_t j = ig[i]; j < ig[i + 1]; j++)
						sum += Matrix[j] * x[jg[j]];
					x[i] = (R[i] - sum) / Diagonal[i];
				}
				return x;
			}
			std::vector<T> Down(std::vector<T> Matrix, std::vector<T> R)
			{
				std::vector<T> x(N);
				for (size_t i = 0; i < N; i++)
				{
					T sum = 0;
					for (size_t j = ig[i]; j < ig[i + 1]; j++)
						sum += Matrix[j] * x[jg[j]];
					x[i] = R[i] - sum;
				}
				return x;
			}

			//Обратный ход
			std::vector<T> Up(std::vector<T> Matrix, std::vector<T> Diagonal, std::vector<T> R)
			{
				for (int i = N - 1; i >= 0; i--)
				{
					R[i] /= Diagonal[i];
					for (size_t j = ig[i]; j < ig[i + 1]; j++)
					{
						size_t p = jg[j];
						R[p] -= Matrix[j] * R[i];
					}
				}
				return R;
			}
			std::vector<T> Up(std::vector<T> Matrix, std::vector<T> R)
			{
				for (int i = N - 1; i >= 0; i--)
				{
					for (size_t j = ig[i]; j < ig[i + 1]; j++)
					{
						size_t p = jg[j];
						R[p] -= Matrix[j] * R[i];
					}
				}
				return R;
			}

			//Умноженик числа на вектор
			std::vector<T> MultVectorOnT(std::vector<T> a, T b)
			{
				for (size_t i = 0; i < a.size(); i++)
					a[i] *= b;
				return a;
			}

			//Mult matrix and vector f
			///<param name = "f"> Matrix be multtiplication on thix vector.</param>
			std::vector<T> MultMatrixVector(std::vector<T> f)
			{
				std::vector<T> v(N);
				for (size_t i = 0; i < N; i++)
				{
					for (size_t j = ig[i]; j < ig[i + 1]; j++)
					{
						v[i] += ggl[j] * f[jg[j]];
						v[jg[j]] += ggu[j] * f[i];
					}
					v[i] += di[i] * f[i];
				}
				return v;
			}
			std::vector<T> MultTMatrixVector(std::vector<T> f)
			{
				std::vector<T> v(N);
				for (size_t i = 0; i < N; i++)
				{
					for (size_t j = ig[i]; j < ig[i + 1]; j++)
					{
						v[i] += ggu[j] * f[jg[j]];
						v[jg[j]] += ggl[j] * f[i];
					}
					v[i] += di[i] * f[i];
				}
				return v;
			}

			//Residual a - b
			std::vector<T> Residual(std::vector<T> a, std::vector<T> b)
			{
				std::vector<T> v;
				v.resize(a.size());
				for (size_t i = 0; i < a.size(); i++)
					v[i] = a[i] - b[i];
				return v;
			}
			//Sum a + b
			std::vector<T> Sum(std::vector<T> a, std::vector<T> b)
			{
				std::vector<T> v;
				v.resize(a.size());
				for (size_t i = 0; i < a.size(); i++)
				{
					v[i] = a[i] + b[i];
				}
				return v;
			}

			//Норма вектора а
			inline T Norm(std::vector<T> a)
			{
				return sqrt(Scolar(a));
			}
			//Скалярное произведение 1 элемента
			inline T Scolar(std::vector<T> a)
			{
				T sum = 0;
				for (size_t i = 0; i < N; i++)
					sum += pow(a[i], 2);
				return sum;
			}
			//Скалярное произведение 2 элемента
			inline T Scolar(std::vector<T> a, std::vector<T> b)
			{
				T sum = 0;
				for (size_t i = 0; i < N; i++)
					sum += a[i] * b[i];
				return sum;
			}

			//factorizationLU is a method in the Matrix class. МБ даже работает
			///<param name = "h"> if this param true, then LU, else diagonal LU </param>
			void factorizationLU(bool h)
			{
				if (h)
				{
					L.resize(ggl.size());
					D.resize(di.size());
					U.resize(ggu.size());
					for (size_t i = 0; i < N; i++)
					{
						T sum = 0;
						for (size_t j = ig[i]; j < ig[i + 1]; j++)
						{
							T sum1 = 0, sum2 = 0;
							int jj = jg[j];
							for (size_t k = ig[i], k2 = ig[jj]; k < j && k2 < ig[jj + 1];)
							{
								int p1 = jg[k];
								int p2 = jg[k2];
								if (p1 == p2)
								{
									sum1 += L[k] * U[k2];
									sum2 += L[k2] * U[k];
									k++; k2++;
								}
								else
								{
									if (p1 < p2)
										k++;
									else
										k2++;
								}
							}
							L[j] = (ggl[j] - sum1) / D[jg[j]];
							U[j] = ggu[j] - sum2;
							sum += L[j] * U[j];
						}
						D[i] = di[i] - sum;
					}
				}
				else
				{
					L.resize(ggl.size());
					U.resize(ggu.size());
					D = di;
				}
			}
			//factorizationLU is a method in the Matrix class. МБ даже работает
			///<param name = "h"> if this param true, then Holisskigo, else diagonal Holisskigo </param>
			void factorizationLLT()
			{
				L.resize(ggl.size());
				D.resize(N);
				for (size_t i = 0; i < N; i++)
					D[i] = sqrt(di[i]);
			}
		};
	}
}


