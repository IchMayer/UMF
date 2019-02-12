#pragma once
#include"Includes.h"
#include "matrix (2).cpp"

using namespace std;

template<class T>
class Matrix
{
public:
	Matrix() { maxiter = 1000;  d = k = n = m = Error = 0; e = 1e-14; w = 1; }
	~Matrix(){}

	//Матрица A, хранящаяся как одномерный массив
	T matrix(size_t i,size_t j);
	//х Уравнения Ax = F
	vector<T> result;
	//Левая часть равенство Ах = F, (F - f) - нвязка
	vector<T> F, f;
	//Инжексы диагоналей (общий случий)
	vector<int> Index;
	size_t k;
	//эпсиол и омега
	T e, w;
	//Максимальное количество итерация
	size_t maxiter;

	//Получение информации об ошибке
	size_t GetLastError();

	void GetResult();

	void SetMatrix(size_t n, size_t last, vector<vector<double>> L, vector<double> f);

	//Открытие файлов
	void OpenMatrix();
	//Сохранение результата
	void SaveResult();

	//Метод Якоби
	void JacobiMethod();
	//Метода Зиделя
	void GaussSeidelMethod();
	//Метод релаксации с квадратным блоком стороной NBlock
	void BlockRelax(int Nblock);
private:
	//Факторизация
	vector<T> Fractorization(int Nblock);
	//Норма вектора (L2)
	inline T NormVector(vector<T>);
	//Умножение матрицы matrix на вектор X
	vector<T> MultMatrixAandVector(vector<T> X);
	//Разность векторов a и b
	vector<T> ResidualVectors(vector<T> a, vector<T> b);
	//Добавление результата для метода Зиделя
	void Gaussnewvector(vector<T> &f, T x, int i);
	//Добавление результата для метода релаксации
	void Blocknewvector(vector<T> &f, T x, int i, int Max);
	//Проверка на выход
	inline bool CheckEnd(); 
	//Функция параметра релаксивности
	T RelaxParam(T x2, T x1);

#ifdef _DEBUG
	T a[100000]; //Максимальный размер 1000 на 100
	const int MAX = 100000;
#else
	T a[100000]; //Максимальный размер 1000 на 100
	const int MAX = 100000;
#endif // DEBUG
	vector<T> buf;
	size_t n, m;
	size_t Error;
	size_t d;
};
/*Errors:
  0 - not error
  1 - matrix overflow
  2 - end method
*/
