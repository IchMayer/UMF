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
	ConvertShape2Index();
	CreateMatrixL();
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

void Model::ConvertShape2Index()
{
	if (!n && !m)
		return;

	index.maxX = m;
	index.maxY = n;

	index.YT = (int)((shape.YT - shape.minY) / (shape.maxY - shape.minY) * n + 0.5);  //0.5 к ближайшему 
	
	// сделать разный шаг

	index.iXT = (int)((shape.iXT - shape.minX) / (shape.maxX - shape.minX) * n + 0.5);  //0.5 к ближайшему
	index.aXT = (int)((shape.aXT - shape.minX) / (shape.maxX - shape.minX) * n + 0.5);  //0.5 к ближайшему

	//сделать разный шаг

	hHI = (shape.maxY - shape.minY) / n;
	hHJ = (shape.maxX - shape.minX) / m;

	//Можно в формулы подставить hHI и hHj будет немного быстрее, но непонятно откуда формулы
}

void Model::CreateMatrixL()
{
	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < m; j++)
			f.push_back(F(i * hHI, j * hHJ));

	for (size_t i = 0; i < n * m; i++)
	{
		//if(Внктри фигуры)
		//if(Вне фигуры)
		//if(Граница)
		//	shitch(какая граница 1 2 или 3)
	}
}
