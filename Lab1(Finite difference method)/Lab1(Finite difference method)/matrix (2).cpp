#include"Includes.h"
#include "matrix (2).h"

template<class T>
matrixLDU<T>::matrixLDU()
{
	Error = NULL;
	n = 0;
	m = 0;
	if (is_same<T, float>::value)
	{
		CompareNumber = pow(10, 6);
		return;
	}
	if (is_same<T, double>::value)
	{
		CompareNumber = pow(10, 14);
		return;
	}
	CompareNumber = 0;
}

template<class T>
matrixLDU<T>::~matrixLDU()
{
}

//Постройка матриц LDU из матрицы A( Ax = F -> LDUx = F)
template<class T>
void matrixLDU<T>::CountLDU(vector<T> &L, vector<T> &D, vector<T> &U, vector<size_t> &ia, TypeMatrix Type, size_t &error)
{
	if (error)
		return;
	switch (Type)
	{
	case TypeMatrix::Prof:
	{
		Error = &error;
		n = D.size();
		m = L.size();
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
			if (abs((D[i] - sum)) < abs(D[i]/CompareNumber))
			{
				*Error = MatrixNotCount;
				return;
			}
			D[i] -= sum;
		}
		break;
	}
	default:
		break;
	}
}

//Поиск X зная LDU 
template<class T>
void matrixLDU<T>::CountX(vector<T> &L, vector<T> &D, vector<T> &U, vector<T> &F, vector<T> &X, vector<size_t> &ia)
{
	if (!Error || *Error)
		return;
	//Поиск Y
	vector<T> &Y = X;
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
	vector<T> &Z = X;
	for (size_t i = 0; i < n; i++)
		Z[i] /= D[i];

	//Поиск X
	for (int i = n - 1; i >= 0; i--)
		for (int j = ia[i + 1] - ia[i] - 1,
			j1 = i - j - 1,
			j2 = ia[i + 1] - j - 1;
			j >= 0 ; j--)
			X[j1++] -= X[i] * U[j2++];
}

template<class T>
void matrixLDU<T>::SaveFile(vector<T> &X)
{
	if (!Error || *Error)
		return;
	ofstream oX("X.txt");
	if (!oX)
	{
		*Error = ErrorOutputFile;
		return;
	}
	for (size_t i = 0; i < X.size(); i++)
		oX << X[i] << endl;
}

template<class T>
void matrixLDU<T>::SaveLDU(vector<T> &L, vector<T> &D, vector<T> &U, TypeMatrix Type)
{
	if (!Error || *Error)
		return;
	switch (Type)
	{
	case TypeMatrix::Prof:
	{
		ofstream oL("L.txt");
		for (size_t i = 0; i < L.size(); i++)
			oL << L[i] << endl;
		oL.close();

		ofstream oU("U.txt");
		for (size_t i = 0; i < U.size(); i++)
			oU << U[i] << endl;
		oU.close();

		ofstream oD("D.txt");
		for (size_t i = 0; i < D.size(); i++)
			oD << D[i] << endl;
		oD.close();

		break;
	}
	default:
		break;
	}
	return 1;
}

template<class T>
matrixA<T>::matrixA()
{
	if (is_same<T, float>::value)
	{
		CompareNumber = pow(10, 6);
		return;
	}
	if (is_same<T, double>::value)
	{
		CompareNumber = pow(10, 14);
		return;
	}
	CompareNumber = 0;
}

template<class T>
matrixA<T>::~matrixA()
{
}


/*
	HighElemGausStraight - Функция поиска решения СЛАУ методом Гаусса с выбором главного коэфициента.
	Матрица будет приводиться к верхнему треугольному виду
*/
template<class T>
void matrixA<T>::CountUpTrianMatrix()
{
	if (Error)
		return;
	for (size_t i = 0; i < n; i++)
	{
		//Опредиление главного Элемента в столбце
		T Max = A[i][i];
		size_t index = i;
		for (size_t j = i+1; j < n; j++)
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
			Error = MatrixNotCount;
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
			if(!A[j][i])
				continue;
			T c = A[j][i]/Max;
			A[j][i] = 0.;
			for (size_t l = i + 1; l < n; l++)
				A[j][l] -= A[i][l] * c;
			F[j] -= F[i] * c;
		}
	}
}

template<class T>
void matrixA<T>::GausBack()
{
	if (Error)
		return;
	for (int i = n - 1; i >= 0; i--)
	{
		double sum = 0;
		for (size_t j = i + 1; j < n; j++)
			sum += A[i][j] * X[j];
		X[i] = (F[i] - sum)/A[i][i];
	}
}

template<class T>
void matrixA<T>::Prof2Full()
{
	if (Error)
		return;
	if (!A.size())
		if (!di.size())
		{
			Error = MatrixNotCount;
			return;
		}
		else
		{
			A.resize(n);
			for (size_t i = 0; i < n; A[i++].resize(n));
		}
		
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

template<class T>
void matrixA<T>::OpenFiles(TypeMatrix typefile, char* path)
{
	switch (typefile)
	{
	case TypeMatrix::Full:
	{
		ifstream ifs("input.txt");
		if (!ifs) {
			Error = ErrorInputFile;
			return;
		}
		size_t n;
		ifs >> n;

		A.reserve(n);

		for (size_t i = 0; i < n; i++)
		{
			A[i].reserve(n);
			for (size_t j = 0; j < n; j++)
				ifs >> A[i][j];
		}

		ifs.close();
		break;
	}
	case TypeMatrix::Prof:
	{
		size_t m;
		vector<size_t> info(2);
		OpenFile("info.txt", path, info, 2);
		n = info[0], m = info[1];
		X.resize(n);
		OpenFile("al.txt", path, al, m);
		OpenFile("au.txt", path, au, m);
		OpenFile("di.txt", path, di, n);
		OpenFile("ia.txt", path, ia, n + 1);
		OpenFile("F.txt", path, F, n);
		break;
	}
	default:
		break;
	}
}

template<class T>
void matrixA<T>::OpenFiles(TypeMatrix typefile, char *path, int k)
{
	switch (typefile)
	{
	case TypeMatrix::Full:
	{
		ifstream ifs("input.txt");
		if (!ifs) {
			Error = ErrorInputFile;
			return;
		}
		size_t n;
		ifs >> n;

		A.reserve(n);

		for (size_t i = 0; i < n; i++)
		{
			A[i].reserve(n);
			for (size_t j = 0; j < n; j++)
				ifs >> A[i][j];
		}

		ifs.close();
		break;
	}
	case TypeMatrix::Prof:
	{
		size_t m;
		vector<size_t> info(2);
		OpenFile("info.txt", path, info, 2);
		n = info[0], m = info[1];
		X.resize(n);
		OpenFile("al.txt", path, al, m);
		OpenFile("au.txt", path, au, m);
		OpenFile("di.txt", path, di, n);
		OpenFile("ia.txt", path, ia, n + 1);
		OpenFile("F.txt", path, F, n);

		if (Error)
			return;

		double p = pow(10, -k);
		di[0] += p;
		F[0] += p;
		break;
	}
	default:
		break;
	}
	return;
}

template<class T>
template<class U>
bool matrixA<T>::OpenFile(const char name[], char *path, vector<U> &a, size_t len)
{
	char Path[255];
	strcpy_s(Path, path);
	strcat_s(Path, name);
	ifstream ia(Path);
	if (!ia)
	{
		Error = ErrorInputFile;
		return false;
	}
	a.resize(len);
	for (size_t i = 0; i < len; i++)
		ia >> a[i];
	ia.close();
	return true;
}

template<class T>
void matrixA<T>::SaveFile()
{
	if (Error)
		return;
	ofstream oX("X.txt");
	if (!oX)
	{
		Error = ErrorOutputFile;
		return;
	}
	for (size_t i = 0; i < X.size(); i++)
		oX << X[i] << endl;
}