#include "Model.h"
#include "Includes.h"


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


	for (size_t i = 0; i < L.size(); i++)
	{
		for (size_t j = 0; j < 5; j++)
		{
			cout << L[i][j] << " ";
		}
		cout << endl;
	}
	int o;
	cin >> o;
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
				f.push_back(0);
			else
			{
				if (i == index.maxY || !i || !j || j == index.maxX || shape.type == shape.T && (index.YT == i || index.iXT == j || index.aXT == j))
				{
					int k;
					
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
						if (!j && (i != index.maxY && i != index.YT ||
							borders[0].type == borders[0].First ||
							i == index.maxY && borders[3].type != borders[3].First ||
							!i && borders[4].type != borders[4].First))
							k = 0;
						if (!i)
							k = 1;
						if (j == index.maxX)
							k = 2;
						if (i == index.maxX)
							k = 3;
						if (i == index.YT && j < index.iXT || i == index.YT)
							k = 4;
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
