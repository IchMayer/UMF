#pragma once
#include"Includes.h"

using namespace std;

class Model
{
public:
	//Format matrix 
	// T - it's 
	// *  *  *  * 
	// *  *  *  *    
	//    *  *
	//    *  *
	// or more point
	enum Format
	{
		Full,
		T
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

		Type type;			
		function<double(double, double)> borderF;
	};

	Model();
	~Model();

	//Create Model

	/// <param name="N">line</param>
	/// <param name="M">column</param>
	/// <param name="format">format of grid</param>
	void Init(size_t N = 4, size_t M = 4, Format format = T);

	/// <param name="border">border</param>
	void AddBorder(Border border);


	
private:
	size_t n, m;									//Matrix n * m
	vector<vector<double>> U;						//U function
	vector<function<double(double, double)>> F;		//F function
	vector<double> f;								//f - local F
	vector<double> u;								//u - local U
	vector<vector<double>> L;						//L 5 diagonal matrix Lu = f
	double hI, hJ;									//Step
	vector<Border> borders;							//Border 
};

