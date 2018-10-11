#include "Solvers.h"

CSolvers::CSolvers(int num_cells_x, int num_cells_y)
{
	R = new double*[num_cells_x];
	U = new double*[num_cells_x];
	V = new double*[num_cells_x];
	P = new double*[num_cells_x];

	for(size_t i = 0; i < num_cells_x; i++)
	{
		R[i] = new double*[num_cells_y];
		U[i] = new double*[num_cells_y];
		V[i] = new double*[num_cells_y];
		P[i] = new double*[num_cells_y];
	}
}

CSolvers::~CSolvers(void)
{
	for(size_t i = 0; i < num_cells_x; i++)
	{
		delete []R[i];
		delete []U[i];
		delete []V[i];
		delete []P[i];
	}

	delete []R;
	delete []U;
	delete []V;
	delete []P;
}

void CSolvers::linear_solver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss)
{
	double *C = new double[numcells];
	double *RC = new double[numcells];
	double *H = new double[numcells];

#pragma omp parallel for simd /*schedule(simd:static)*/ schedule(dynamic,64) 
		for (int i = 1; i < numcells; i++)
		{

			C[i - 1] = sqrt(GAMMA*P[i - 1] / R[i - 1]);
			C[i] = sqrt(GAMMA*P[i] / R[i]);

			RC[i - 1] = R[i - 1] * C[i - 1];
			RC[i] = R[i] * C[i];

			H[i - 1] = 1.0 / (RC[i - 1]);
			H[i] = 1.0 / (RC[i]);

			if (U[i - 1] > C[i - 1])
			{
				pss[i] = P[i - 1];
				uss[i] = U[i - 1];
				dss[i] = R[i - 1];

			}
			else if (U[i] < -C[i])
			{
				pss[i] = P[i];
				uss[i] = U[i];
				dss[i] = R[i];
			}
			else
			{
				pss[i] = (U[i - 1] - U[i] + P[i - 1] * H[i - 1] + P[i] * H[i]) / (H[i - 1] + H[i]);
				uss[i] = (RC[i - 1] * U[i - 1] + RC[i] * U[i] + P[i - 1] - P[i]) / (RC[i - 1] + RC[i]);

			    //if (uss[i] > 0) dss[i] = R[i - 1] - R[i - 1] / C[i - 1] * (uss[i] - U[i - 1]);
			    //else dss[i] = R[i] + R[i] / C[i] * (uss[i] - U[i]);

				if (uss[i] > 0) dss[i] = R[i - 1] * (1 + (pss[i] - P[i - 1])/(RC[i - 1] * C[i - 1]));
			    else dss[i] = R[i] * (1 + (pss[i] - P[i])/(RC[i] * C[i]));
			}

		}

	delete[] RC;
	delete[] C;
	delete[] H;
}

void CSolvers::solve(double*** V_init, CSolvers* mesh)
{
	double num_cells[2];

	num_cells[0] = mesh->getNumCells()[0];
	num_cells[1] = mesh->getNumCells()[1];

	R = V_init[0];
	U = V_init[1];
	V = V_init[2];
	P = V_init[3];

	linearSolver(num_cells, R, U, P, dss, uss, pss);

	double **UFLUX, **VFLUX, **FR, **FRU, **FRV, **FRE;

#pragma omp for simd schedule(simd:static) 
	for (int i = 0; i <= numcells; i++)
	{
		UFLUX[i] = uss[i];
		FR[i] = dss[i] * uss[i];
		FRU[i] = dss[i] * uss[i] * uss[i] + pss[i];
		FRE[i] = (pss[i] / (GAMMA - 1.0) + 0.5*dss[i] * uss[i] * uss[i])*uss[i] + pss[i] * uss[i];
	}
}