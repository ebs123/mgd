// godunov.cpp : Defines the entry point for the console application.
//

#include "includes\Solvers.h"
//fucking govnocode)))))
#include <fstream>
#include <vector>
#include "includes\flow_parameters.h"

const std::vector<std::string> explode(const std::string& s, const char& c)
{
	std::string buff("");
	std::vector<std::string> v;
	
	for(size_t i = 0; i < s.length(); i++)
	{
		if(s[i] != c) buff+=s[i]; else
		if(s[i] == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if(buff != "") v.push_back(buff);
	
	return v;
}

int main(int argc, char* argv)
{
	CSolvers *solver = new CSolvers;
	int n_steps = 200000;
	int num_cells[2];
	int n_save = 1000000;
	double ***V_init;
	int problem_dimension = 2;
	double *domain_length = new double[problem_dimension];
	int *mesh_size = new int[problem_dimension];
	bool continue_calc = false;

	domain_length[0] = 2 * pi;
	domain_length[1] = 2 * pi;
	mesh_size[0] = 250.;
	mesh_size[1] = 250.;

	COutput *output = new COutput;
	CMeshGenerator *mesh = new CMeshGenerator(domain_length, mesh_size, problem_dimension);
	CFlowParameters *flow_parameters = new CFlowParameters(mesh);

	num_cells[0] = mesh->getNumCells()[0];
	num_cells[1] = mesh->getNumCells()[1];

	//NTracer::openFile("flowTrace.dat");
	V_init = new double**[4];//r,u,v,p
	
	for(size_t i = 0; i < 4; i++)
	{
		V_init[i] = new double*[num_cells[0]];
	}
	for(size_t i = 0; i < 4; i++)
		for(size_t j = 0; j < num_cells[0]; j++)
		{
			V_init[i][j] = new double[num_cells[1]];
		}

	double *x = mesh->getMeshComponent(0);
	double *y = mesh->getMeshComponent(1);

	if(continue_calc)
	{
		std::ifstream file_stream(PATH_TO_DATA_FILE);
		std::string file_data_line;
		int i = 0;
		int j = 0;

	    if(file_stream.is_open())
		{
			while (getline(file_stream, file_data_line))
			{
				std::vector<std::string> data_arr = explode(file_data_line, ' ');
				if(data_arr.size() == 7)
				{
					V_init[0][i][j] = atof(data_arr[2].c_str());
					V_init[1][i][j] = atof(data_arr[3].c_str());
					V_init[2][i][j] = atof(data_arr[4].c_str());
					V_init[3][i][j] = atof(data_arr[5].c_str());
					j++;
					if(j == mesh->getNumCells()[1])
					{
						j = 0;
						i++;
					}

					if(i == mesh->getNumCells()[0])
					{
						break;
					}
				}
			}
		}
		else
		{
			std::cout << "Can't open file " << PATH_TO_DATA_FILE << " for read";
		}
		file_stream.close();
	}
	else
	{
		for(size_t i = 0; i < num_cells[0]; i++)
			for(size_t j = 0; j < num_cells[1]; j++)
			{
				V_init[0][i][j] = 10;
				V_init[1][i][j] = .1 * sin(2 * y[j]);// + .01 * sin(x[i]);
				V_init[2][i][j] = - .1 * sin(2 * x[i]);// + .01 * sin(y[j]);
				V_init[3][i][j] = 10000;
				/*if(y[j] < 0.5 + .01 * cos(6 * pi * x[i]))
				{
					V_init[0][i][j] = 1.;
					V_init[1][i][j] = 0;
					V_init[2][i][j] = 0;
					V_init[3][i][j] = .1 * (1.5 - y[j]) + .01;
				}
				else
				{
					V_init[0][i][j] = 2.;
					V_init[1][i][j] = 0;
					V_init[2][i][j] = 0;
					V_init[3][i][j] = 2. * .1 * (1 - y[j]) + .01;
				}*/
				//if(x[i] <= .1)
				//{
				//	V_init[0][i][j] = 1.2714;
				//	V_init[1][i][j] = 0.2928;
				//	V_init[2][i][j] = 0.;
				//	V_init[3][i][j] = 1.4017;
				//}
				//else
				//{
				//	V_init[0][i][j] = 1.;
				//	V_init[1][i][j] = 0.;
				//	V_init[2][i][j] = 0.;
				//	V_init[3][i][j] = 1.;
				//}
			}
	}

	output->save2dPlot("initial.dat", mesh, 0, V_init[0], V_init[1], V_init[2], V_init[3], flow_parameters->getFlowEntropy(V_init[0], V_init[3]));
	//output->save1dPlotXAxis("initialX.dat", mesh, 0, 50, V_init[0], V_init[1], V_init[2], V_init[3]);
	//output->save1dPlotYAxis("initialY.dat", mesh, 0, 50, V_init[0], V_init[1], V_init[2], V_init[3]);
	//exit(1);

	solver->solve(V_init, n_steps, n_save, mesh, SLIP, SLIP);
	return 0;
}