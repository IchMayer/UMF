#pragma once

#include"Includes.h"

using namespace std;

enum class TypeMatrix
{
	Full,
	Prof,
	Tape
};

enum {
	NoErrors = 0,
	ErrorInputFile = 1,
	ErrorOutputFile = 2,
	MatrixNotCount = 3
};

template<class T>
class matrixA
{
public:

	vector<vector<T>> A;	//Матрица в полном формате
	vector<T> al, au, di;	//Матрица в профильном формате
	vector<size_t> ia;		//Матрица в профильном формате

	vector<T> F;			//Мектор правой стороны СЛАУ Ax = F
	vector<T> X;			//Результат СЛАУ 
	size_t Error;			//Информация об ошибках

	matrixA();
	~matrixA();				

	void SetN(int N) { n = N; }
	void SaveFile();				//Сохранение решения СЛАУ в текстоывый документ
	void CountUpTrianMatrix();		//Постройка верхней треугольной матрицы
	void GausBack();				//Обратный метод Гаусса
	void Prof2Full();			    //Перевод из профильного формата в полный

	void OpenFiles(TypeMatrix, char *);    //Открытие файлов
	void OpenFiles(TypeMatrix, char *, int k); 
private:
	double CompareNumber;
	size_t n; //Размерность матрицы

	template <class U>
	bool OpenFile(const char name[], char *path, vector<U> &a, size_t len);
};

template<class T>
class matrixLDU
{
public:
	matrixLDU();
	~matrixLDU();

	void SetNErrorM(size_t N, size_t M, size_t &error) { n = N; m = M; Error = &error; }
	void CountLDU(vector<T> &L, vector<T> &D, vector<T> &U, vector<size_t> &ia, TypeMatrix type, size_t &Error);//Перевод из профильного формата matrixA в профильный формат LDU
	void CountX(vector<T> &L, vector<T> &D, vector<T> &U, vector<T> &F, vector<T> &X, vector<size_t> & ia);		//Поиск решения уравнения LDUx = F
	void SaveFile(vector<T> &);																	//Сохранения результата СЛАУ LDU x = F в файл
	void SaveLDU(vector<T> &L, vector<T> &D, vector<T> &U, TypeMatrix Type);						//Сохранения LDU в файл
private:
	double CompareNumber;
	size_t n; //Размерность матрицы
	size_t m; //Количество недиагональных элементов в матрице
	size_t *Error; //Информация об ошибках
};