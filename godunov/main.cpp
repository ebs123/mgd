// godunov.cpp : Defines the entry point for the console application.
//

#include "includes\Solvers.h"


int main(int argc, char* argv)
{
	CSolvers *solver = new CSolvers;
	int n_steps = 10;
	int num_cells[2];
	int n_save = 1;
	double ***V_init;
	int problem_dimension = 2;
	double *domain_length = new double[problem_dimension];
	int *mesh_size = new int[problem_dimension];

	domain_length[0] = 1.;
	domain_length[1] = 1.;
	mesh_size[0] = 100.;
	mesh_size[1] = 100.;

	COutput *output = new COutput;
	CMeshGenerator *mesh = new CMeshGenerator(domain_length, mesh_size, problem_dimension);
	num_cells[0] = mesh->getNumCells()[0];
	num_cells[1] = mesh->getNumCells()[1];

	NTracer::openFile("flowTrace.dat");
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

	for(size_t i = 0; i < num_cells[0]; i++)
		for(size_t j = 0; j < num_cells[1]; j++)
		{
			V_init[0][i][j] = 1.;
			V_init[1][i][j] = 0.;
			V_init[2][i][j] = 0.;
			if(x[i] < domain_length[0] * .5)
			{
				V_init[3][i][j] = 1.;
			}
			else
			{
				V_init[3][i][j] = 2.;
			}
		}

	output->save2dPlot("initial.dat", mesh, V_init[0], V_init[1], V_init[2], V_init[3], V_init[3]);

	solver->solve(V_init, n_steps, n_save, mesh);
	return 0;
}
