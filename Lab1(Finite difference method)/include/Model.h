#pragma once
#include"Includes.h"

using namespace std;

class Model
{
public:
	//Форма
	struct Shape
	{
		//Format matrix 
		// T - it's 
		// *  *  *  *	<- maxY
		// *  *  *  *   <- YT
		//    *  *  
		//    *  *		<- minY
		// A  A  A  A
		// |  |  |  |
		// m  i  a  m
		// i  X  X  a
		// n  T  T  x
		// X        X
		// or more point
		enum Format
		{
			Rectangle,
			T
		};
		Format type;
		double minX, maxX;
		double minY, maxY;
		double iXT, aXT;
		double YT;
	};

	//Граница
	struct Border
	{
		enum Type
		{
			First,
			Second,
			Third
		};

		double b;
		Type type;			
		function<double(double, double)> borderF;
	};

	Model();
	~Model();

	//Create Model

	/// <param name="N">line</param>
	/// <param name="M">column</param>
	/// <param name="format">format of grid</param>
	void CreateGrid(size_t N = 4, size_t M = 4);

	/// <param name="shape">shape of area</param>
	void SetShape(Shape shape);

	/// <param name="border">border</param>
	void AddBorder(Border border);

	/// <param name="F">F from Lu = F. </param>
	void SetF(function<double(double, double)> F);

	/// <param name="lamda">New lamda </param>
	void SetLamda(double lamda);

	/// <param name="gamma">New gamma </param>
	void SetGamma(double gamma);


	void FDM();
	
	vector<vector<double>> GetResult();

	/// <param name="Kx">Коэфициент неравномерности по x </param>
	/// <param name="Ky">Коэфициент неравномерности по y </param>
	void SetKGrid(double Kx, double Ky);

	double GethHI();
	double GethHJ();

private:
	struct ShapeIndex
	{
		size_t maxX;
		size_t maxY;
		size_t iXT, aXT;
		size_t YT;
	};

	double lamda;									//lamda
	double gamma;									//gamma
	size_t n, m;									//Matrix n * m
	function<double(double, double)> F;				//F function
	vector<double> f;								//f - local F
	vector<double> u;								//u - local U
	vector<vector<double>> L;						//L 5 diagonal matrix Lu = f
	double hHI, hHJ;								//первый шаг сетки
	double hLT, hJT;								//коэффициент неравномерности 
	Shape shape;									//Shape of area
	ShapeIndex index;								//Shape in index
	vector<Border> borders;							//Border of area

	void ConvertShape2Index();		
	void ConvertNeravShape2Index();

	void ConvertL2MatrixFormat();					//Matrix - программа из прошлого года, которая решает диагональные матрицы методлом якоби и зегеля

	vector<vector<double>> ConvertL2FullMatrix();

	void CreateMatrixL();
};

