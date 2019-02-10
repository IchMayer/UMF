#include "Model.h"
#include "Includes.h"


Model::Model()
{
}


Model::~Model()
{
}

void Model::CreateGrid(size_t N, size_t M)
{
	n = N;
	m = M;
}

void Model::SetShape(Shape s)
{
	shape = s;
}

void Model::AddBorder(Border border)
{
	borders.push_back(border);
}

void Model::SetF(function<double(double, double)> f)
{
	F = f;
}

void Model::FDM()
{
}

vector<vector<double>> Model::GetResult()
{	
	//Conver line to matrix
	vector<vector<double>> U;
	for (size_t i = 0; i < n; i++)
	{
		for (size_t i = 0; i < m; i++)
		{

		}
	}
	return U;
}
