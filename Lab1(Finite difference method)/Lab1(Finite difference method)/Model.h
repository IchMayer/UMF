#pragma once
#include"Includes.h"

using namespace std;

class Model
{
public:
	enum Format
	{
		Full,
		T
	};
	Model();
	~Model();

	//Create Model
	void init();

	/// <param name="n">column</param>
	/// <param name="m">line</param>
	/// <param name="format">format of grid</param>
	void init(size_t n, size_t m, Format format);


private:
	size_t n, m;					// Matrix n * m
	vector<vector<double>> U;		//U function
	double hI, hJ;					//Step
};

