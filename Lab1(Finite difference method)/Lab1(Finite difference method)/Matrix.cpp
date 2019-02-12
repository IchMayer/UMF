#include"Includes.h"
#include"Matrix.h"

template<class T>
inline T Matrix<T>::matrix(size_t i, size_t j)
{
	if (!n)
		return false;
	return a[i * n + j];
}

template<class T>
size_t Matrix<T>::GetLastError()
{
	return Error;
}

template<class T>
void Matrix<T>::SetMatrix(size_t N, size_t M, vector<vector<double>> L, vector<double> f)
{
	n = N;
	m = 5;

	k = 0;
	Error = 0;

	if (n * m > MAX)
	{
		Error = 1;
		return;
	}


	Index.resize(m);
	F.resize(n);
	result.resize(n, 0);

	Index[0] = -1 * M;
	Index[1] = -1;
	Index[2] = 0;
	Index[3] = 1;
	Index[4] = M;
	d = 2;

	for (size_t i = 0, k = m - 1; i < m; i++, k--)
		for (size_t j = 0, i1 = i * n + j; j < n; j++)
			a[i1++] = L[j][k];

	for (size_t i = 0; i < n; i++)
		result[i] = 0;

	for (size_t i = 0; i < n; i++)
		F[i] = f[i];
}

template<class T>
void Matrix<T>::OpenMatrix()
{
	ifstream index("index.txt");
	index >> n;
	index >> m;

	k = 0;
	Error = 0;

	if (n * m > MAX)
	{
		Error = 1;
		return;
	}

	Index.resize(m);
	F.resize(n);
	result.resize(n, 0);

	for (size_t i = 0; i < m; i++)
	{
		index >> Index[i];
		if (!Index[i])
			d = i;
	}
	index.close();

	ifstream matrix("matrix.txt");
	for (size_t i = 0; i < m; i++)
		for (size_t j = 0, i1 = i * n + j; j < n; j++)
			matrix >> a[i1++];
	matrix.close();

	ifstream X("X.txt");
	if (X)
		for (size_t i = 0; i < n; i++)
			X >> result[i];
	else
		for (size_t i = 0; i < n; i++)
			result[i] = 0;
	X.close();

	ifstream f("F.txt");
	for (size_t i = 0; i < n; i++)
		f >> F[i];
	f.close();
}

template<class T>
void Matrix<T>::SaveResult()
{
	ofstream X("Result.txt");
	for (size_t i = 0; i < n; ++i)
		X << result[i];
}

template<class T>
void Matrix<T>::JacobiMethod()
{
	if (CheckEnd())
	{
		Error = 2;
		return;
	}
	for (size_t i = 0; i < n; i++)
	{
#ifdef _DEBUG
		T i1 = RelaxParam(result[i] + (F[i] - f[i]) / matrix(d, i), result[i]), i2 = matrix(d, i), i3 = F[i] - f[i], i4 = result[i] + (F[i] - f[i]) / matrix(d, i), i5 = result[i];
#endif
		result[i] = RelaxParam(result[i] + (F[i] - f[i]) / matrix(d, i), result[i]);
	}
	k++;
}

template<class T>
void Matrix<T>::GaussSeidelMethod()
{
	if (CheckEnd())
	{
		Error = 2;
		return;
	}
	for (size_t i = 0; i < n; i++)
	{
#ifdef _DEBUG
		T i1 = RelaxParam(result[i] + (F[i] - f[i]) / matrix(d, i), result[i]), i2 = matrix(d, i), i3 = F[i] - f[i];
#endif
		T res = RelaxParam(result[i] + (F[i] - f[i]) / matrix(d, i), result[i]);
		Gaussnewvector(f, res, i);
		result[i] = res;
	}
	k++;
}


//Only 3 deogonal * 3  format
template<class T>
void Matrix<T>::BlockRelax(int Nblock)
{
	vector<matrixA<double>> Fractoriz(n / Nblock + (n % Nblock) / n);
	for (size_t i = 0; i < Fractoriz.size(); i++)
	{
		matrixA<double> &b = Fractoriz[i];
		b.ia.push_back(0);
		for (size_t j = 0, j2 = Nblock * i + j; j < n && j < Nblock; j++, j2++)
		{
			b.di.push_back(matrix(d, j2));
			b.ia.push_back(j);
		}
		b.Error = 0;
		b.SetN(2);
		for (size_t j = 0, d1 = d - 1, j2 = Nblock * i + j; j < n && j < Nblock - 1; j++, j2++)
			b.al.push_back(matrix(d1, j2));
		for (size_t j = 1, d2 = d + 1, j2 = Nblock * i + j; j < n && j < Nblock; j++, j2++)
			b.au.push_back(matrix(d2, j2));
		matrixLDU<double> ldu;
		ldu.CountLDU(b.al, b.di, b.au, b.ia, TypeMatrix::Prof, b.Error);
	}
	while (1)
	{
		if (CheckEnd())
		{
			Error = 2;
			return;
		}
		for (size_t i = 0, i1 = 0, i2 = 0; i < Fractoriz.size(); i++)
		{
			f = MultMatrixAandVector(result);
			matrixA<double> &b = Fractoriz[i];
			vector<double> wR(0);
			for (size_t k = 0; i1 < n && k < Nblock; k++, i1++)
				wR.push_back((F[i1] - f[i1]) * w);
			matrixLDU<double> ldu;
			ldu.SetNErrorM(b.di.size(), b.au.size(), b.Error);
			vector<double> newResult(wR.size());
			ldu.CountX(b.al, b.di, b.au, wR, newResult, b.ia);
			for (size_t k = 0, k2 = min(n - 1, i2 + Nblock - 1); k < Nblock && i2 < n; k++, i2++)
			{
		//		Blocknewvector(f, newResult[k], i2, k2);
				result[i2] += newResult[k];
			}
		}
		k++;
	}
}

template<class T>
vector<T> Matrix<T>::Fractorization(int Nblock)
{
	vector<T> Fractorization(n * Nblock);
	size_t deistv = n / Nblock + (n % Nblock) / n;
	for (size_t i = 0, k1 = 0; i < deistv; i++, k1 += Nblock)
	{
		for (size_t j = 0; j < Nblock && k1 + j < n; j++)
		{
			matrixA<T> m;
		}
	}
}

template<class T>
inline T Matrix<T>::NormVector(vector<T> v)
{
	/*T sum = 0;
	for (size_t i = 0; i < v.size(); i++)
		sum += abs(v[i]);
	return sum;*/
	/*T sum = v[0];
	for (size_t i = 1; i < v.size(); i++)
		if (sum < v[i])
			sum = v[i];
	return sum;*/
	return sqrt(inner_product(v.begin(), v.end(), v.begin(), 0.0));
}

template<class T>
vector<T> Matrix<T>::MultMatrixAandVector(vector<T> X)
{
	vector<T> F(n);
	for (size_t i = 0; i < m; i++)
		if (Index[i] >= 0)
			for (size_t j = Index[i], j1 = j - Index[i]; j < n; j++)
			{
#ifdef _DEBUG
				T i1 = X[j], i2 = matrix(i, j), i3 = F[j - Index[i]];
#endif
				F[j1++] += X[j] * matrix(i, j);
			}
		else
			for (size_t j = 0, j1 = j - Index[i]; j - Index[i] < n; j++)
			{
#ifdef _DEBUG
				T i1 = X[j], i2 = matrix(i, j), i3 = F[j - Index[i]];
#endif
				F[j1++] += X[j] * matrix(i, j);
			}
	return F;
}

template<class T>
vector<T> Matrix<T>::ResidualVectors(vector<T> a, vector<T> b)
{
	vector<T> c(a.size());
	for (int i = 0; i < n; i++)
		c[i] = a[i] - b[i];
	return c;
}

template<class T>
void Matrix<T>::Gaussnewvector(vector<T>& vect, T x, int j)
{
	T Res = x - result[j];
	for (size_t i = 0; i < m; i++)
	{
		int j1 = j - Index[i];
		if (j1 > j && j1 < n)
		{
#ifdef _DEBUG
			T i1 = vect[j1], i2 = matrix(i, j);
#endif
			vect[j1] += Res * matrix(i, j);
		}
	}
}

template<class T>
void Matrix<T>::Blocknewvector(vector<T>& vect, T x, int j, int MMAX)
{
	T Res = x - result[j];
	for (size_t i = 0; i < m; i++)
	{
		int j1 = j - Index[i];
		if (j1 > MMAX && j1 < n)
		{
#ifdef _DEBUG
			T i1 = vect[j1], i2 = matrix(i, j);
#endif
			vect[j1] += Res * matrix(i, j);
		}
	}
}

template<class T>
inline bool Matrix<T>::CheckEnd()
{
	if (k > maxiter)
		return true;
	f = MultMatrixAandVector(result);
	if (NormVector(ResidualVectors(F, f)) / NormVector(F) < e)
	{
		vector<T> j = MultMatrixAandVector(result), i = ResidualVectors(F, MultMatrixAandVector(result)); T i1 = NormVector(i), i2 = NormVector(F), i3 = i1 / i2;
		T sum = 0, sum2 = 0;
		for (size_t i = 0; i < result.size(); i++)
		{
			sum += pow(i + 1, 2);
			sum2 += pow(result[i] - i - 1, 2);
		}
		T coudA = sqrt(sum2) / sqrt(sum) / i3;
		wcout << L"Coud A = " << coudA << endl;
		return true;
	}
	return false;
}

template<class T>
T Matrix<T>::RelaxParam(T x2, T x1)
{
	return w * x2 + (1. - w)*x1;
}

