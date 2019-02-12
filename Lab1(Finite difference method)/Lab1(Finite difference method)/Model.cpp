#include "Model.h"
#include "Includes.h"
#include "Matrix.cpp"


Model::Model()
{
	lamda = 1;
	n = 4;
	m = 4;
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

void Model::SetLamda(double l)
{
	lamda = l;
}

void Model::FDM()
{
	ConvertShape2Index();
	CreateMatrixL();

	Matrix<double> a;

	ConvertL2MatrixFormat();

	a.SetMatrix(L.size(), n, L, f);

	//a.GaussSeidelMethod();
	//a.JacobiMethod();

	do
	{
		a.GaussSeidelMethod();
	} while (a.GetLastError() != 2);

	cout << endl;


	u = a.result;
}

vector<vector<double>> Model::GetResult()
{	
	//Conver line to matrix
	vector<vector<double>> U(n, vector<double>(m));
	for (int i = n - 1, k = 0; i >= 0; i--)
		for (size_t j = 0; j < m; j++, k++)
			U[i][j] = u[k];
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

	//ћожно в формулы подставить hHI и hHj будет немного быстрее, но непон€тно откуда формулы
}

void Model::ConvertL2MatrixFormat()
{
	double buf;
	double size = L.size();
	for (size_t i = 0, j = size - 1; i < size - n; i++, j--)
	{
		L[i][4] = L[i + n][4];
		L[j][0] = L[j - n][0];
	}
	for (size_t i = 0, j = size - 1; i < n; i++, j--)
	{
		L[j][4] = 0;
		L[i][0] = 0;
	}
	for (size_t i = 0, j = size - 1; i < size - 1; i++, j--)
	{
		L[i][3] = L[i + 1][3];
		L[j][1] = L[j - 1][1];
	}
	L[size - 1][3] = 0;
	L[0][1] = 0;
}

void Model::CreateMatrixL()
{
	m++;
	n++;
	L.resize(n * m);



	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < m; j++)
		{
			L[i * m + j].resize(5);
			if (shape.type == shape.T && i < index.YT && (j < index.iXT || j > index.aXT))
			{
				L[i * m + j][2] = 1;
				f.push_back(0);
			}
			else
			{
				if (i == index.maxY || !i || !j || j == index.maxX || shape.type == shape.T && (index.YT == i && (j <= index.iXT || j >= index.aXT) ||
					i <= index.YT && (index.iXT == j || index.aXT == j )))
				{
					int k = 0;
					
					if (shape.type == shape.Rectangle)
					{
						if (!j && (i != index.maxY && i ||
							borders[0].type == borders[0].First ||
							i == index.maxY && borders[3].type != borders[3].First ||
							!i && borders[1].type != borders[1].First))
							k = 0;
						if (!i && (j && j != index.maxX ||
							borders[1].type == borders[1].First ||
							j == index.maxX && borders[2].type != borders[2].First ||
							!j && borders[0].type != borders[0].First))
							k = 1;
						if (j == index.maxX && (i && i != index.maxY || 
							borders[2].type == borders[2].First ||
							i == index.maxY && borders[3].type != borders[3].First ||
							!i && borders[1].type != borders[1].First))
							k = 2;
						if (i == index.maxX && (j && j != index.maxX ||
							borders[3].type == borders[3].First ||
							j == index.maxX && borders[2].type != borders[2].First ||
							!j && borders[0].type != borders[0].First))
							k = 3;

						f.push_back(borders[k].borderF(i * hHI, j * hHJ));

						if (borders[k].type == borders[k].First)
						{
							L[i * m + j][2] = 1;
							continue;
						}

						bool p = borders[k].type == borders[k].Second;

						switch (k)
						{
						case 0:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][3] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 1:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][4] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 2:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][1] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 3:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][0] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						}
						continue;
					}

					if (shape.type == shape.T)
					{
						if (!j && (i != index.maxY && i ||
							borders[0].type == borders[0].First ||
							i == index.maxY && borders[3].type != borders[3].First ||
							i == index.YT && borders[4].type != borders[4].First))
							k = 0;
						if (!i && (j != index.aXT  && j != index.maxX ||
							borders[1].type == borders[1].First ||
							j == index.aXT && borders[6].type != borders[6].First ||
							j == index.iXT && borders[5].type != borders[5].First))
							k = 1;
						if (j == index.maxX && (i != index.YT && i != index.maxY ||
							borders[2].type == borders[2].First ||
							i == index.maxY && borders[3].type != borders[3].First ||
							i == index.YT && borders[7].type != borders[7].First))
							k = 2;
						if (i == index.maxX && (j && j != index.maxX ||
							borders[3].type == borders[3].First ||
							j == index.maxX && borders[2].type != borders[2].First ||
							!j && borders[0].type != borders[0].First))
							k = 3;
						if (i == index.YT && j <= index.iXT && (j && j != index.iXT ||
							borders[4].type == borders[4].First ||
							!j && borders[0].type != borders[0].First ||
							j == index.iXT && borders[5].type != borders[5].First))
							k = 4;
						if (j == index.iXT && i <= index.YT && (i && i != index.YT ||
							borders[5].type == borders[5].First ||
							!i && borders[1].type != borders[1].First ||
							i == index.YT && borders[4].type != borders[4].First))
							k = 5;
						if (j == index.aXT && i <= index.YT && (i && i != index.YT ||
							borders[6].type == borders[6].First ||
							!i && borders[1].type != borders[1].First ||
							i == index.YT && borders[7].type != borders[7].First))
							k = 6;
						if (i == index.YT && j >= index.aXT && (j != index.maxX && j != index.aXT ||
							borders[7].type == borders[7].First ||
							j == index.maxX && borders[2].type != borders[2].First ||
							j == index.aXT && borders[6].type != borders[6].First))
							k = 7;

						f.push_back(borders[k].borderF(i * hHI, j * hHJ));

						if (borders[k].type == borders[k].First)
						{
							L[i * m + j][2] = 1;
							continue;
						}

						bool p = borders[k].type == borders[k].Second;

						switch (k)
						{
						case 0:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][3] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 1:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][4] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 2:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][1] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 3:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][0] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 4:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][4] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 5:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][3] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 6:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][1] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						case 7:
							if (p)
							{
								L[i * m + j][2] = -lamda / hHI;
								L[i * m + j][4] = lamda / hHI;
							}
							else
							{
								//need more code
							}
							break;
						}
						continue;
					}
				}
				else
				{
					L[i * m + j][0] = -lamda / pow(hHI, 2);
					L[i * m + j][1] = -lamda / pow(hHJ, 2);
					L[i * m + j][2] = 2 * lamda *(1.0 / pow(hHJ, 2) + 1.0 / pow(hHI, 2) );
					L[i * m + j][3] = -lamda / pow(hHJ, 2);
					L[i * m + j][4] = -lamda / pow(hHI, 2);

					f.push_back(F(i * hHI, j * hHJ));
				}
			}
		}
	}
}
