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

void CSolvers::linear_solver(int* numcells, double** R, double** U, double** V, double** P, double*** dss, double*** uss, 
							 double*** vss, double*** pss)
{
	double **C = new double[numcells[0]];
	double **RC = new double[numcells[0]];
	double **H = new double[numcells[0]];

	for(size_t i = 0; i < numcells[0]; i++)
	{
		C[i] = new double[numcells[1]];
		RC[i] = new double[numcells[1]];
		H[i] = new double[numcells[1]];
	}

#pragma omp parallel for collapse(2)
	for(int i = 1; i < numcells[0]; i++)
		for(int j = 1; j < numcells[1]; j++)
		{
			C[i - 1][j] = sqrt(GAMMA*P[i - 1][j] / R[i - 1][j]);
			C[i][j] = sqrt(GAMMA*P[i][j] / R[i][j]);
			C[i][j - 1] = sqrt(GAMMA*P[i][j - 1] / R[i][j - 1]);

			RC[i - 1][j] = R[i - 1][j] * C[i - 1][j];
			RC[i][j] = R[i][j] * C[i][j];
			RC[i][j - 1] = R[i][j - 1] * C[i][j - 1];

			H[i - 1][j] = 1.0 / (RC[i - 1][j]);
			H[i][j] = 1.0 / (RC[i][j]);
			H[i][j - 1] = 1.0 / (RC[i][j - 1]);

			//one direction
			if(U[i - 1][j] > C[i - 1][j])
			{
				pss[0][i][j] = P[i - 1][j];
				uss[0][i][j] = U[i - 1][j];
				dss[0][i][j] = R[i - 1][j];
				vss[0][i][j] = V[i - 1][j];
			}
			else if (U[i][j] < -C[i][j])
			{
				pss[0][i][j] = P[i][j];
				uss[0][i][j] = U[i][j];
				dss[0][i][j] = R[i][j];
				vss[0][i][j] = V[i][j];
			}
			else
			{
				pss[0][i][j] = (U[i - 1][j] - U[i][j] + P[i - 1][j] * H[i - 1][j] + P[i][j] * H[i][j]) / 
					(H[i - 1][j] + H[i][j]);
				uss[0][i][j] = (RC[i - 1][j] * U[i - 1][j] + RC[i][j] * U[i][j] + P[i - 1][j] - P[i][j]) / 
					(RC[i - 1][j] + RC[i][j]);

				//if (uss[i] > 0) dss[i] = R[i - 1] - R[i - 1] / C[i - 1] * (uss[i] - U[i - 1]);
				//else dss[i] = R[i] + R[i] / C[i] * (uss[i] - U[i]);

				if (uss[0][i][j] > 0)
				{
					dss[0][i][j] = R[i - 1][j] * (1 + (pss[0][i][j] - P[i - 1][j])/(RC[i - 1][j] * 
						C[i - 1][j]));
					vss[0][i][j] = V[i - 1][j];
				}
				else
				{
					dss[0][i][j] = R[i][j] * (1 + (pss[0][i][j] - P[i][j])/(RC[i][j] * C[i][j]));
					vss[0][i][j] = V[i][j];
				}
			}

			//two direction
			if(V[i][j - 1] > C[i][j - 1])
			{
				pss[1][i][j] = P[i][j - 1];
				uss[1][i][j] = U[i][j - 1];
				dss[1][i][j] = R[i][j - 1];
				vss[1][i][j] = V[i][j - 1];
			}
			else if (V[i][j] < -C[i][j])
			{
				pss[1][i][j] = P[i][j];
				uss[1][i][j] = U[i][j];
				dss[1][i][j] = R[i][j];
				vss[1][i][j] = V[i][j];
			}
			else
			{
				pss[1][i][j] = (V[i][j - 1] - V[i][j] + P[i][j - 1] * H[i][j - 1] + P[i][j] * H[i][j]) / 
					(H[i][j - 1] + H[i][j]);
				vss[1][i][j] = (RC[i][j - 1] * V[i][j - 1] + RC[i][j] * V[i][j] + P[i][j - 1] - P[i][j]) / 
					(RC[i][j - 1] + RC[i][j]);

				//if (uss[i] > 0) dss[i] = R[i - 1] - R[i - 1] / C[i - 1] * (uss[i] - U[i - 1]);
				//else dss[i] = R[i] + R[i] / C[i] * (uss[i] - U[i]);

				if (vss[1][i][j] > 0)
				{
					dss[1][i][j] = R[i][j - 1] * (1 + (pss[1][i][j] - P[i][j - 1])/(RC[i][j - 1] * 
						C[i][j - 1]));
					uss[1][i][j] = U[i - 1][j];
				}
				else
				{
					dss[1][i][j] = R[i][j] * (1 + (pss[1][i][j] - P[i][j])/(RC[i][j] * C[i][j]));
					uss[1][i][j] = U[i][j];
				}
			}

		}

	for(size_t i = 0; i < numcells[0]; i++)
	{
		delete[] RC[i];
		delete[] C[i];
		delete[] H[i];
	}

	delete[] RC;
	delete[] C;
	delete[] H;

}

void CSolvers::solve(double*** V_init, CSolvers* mesh)
{
	int num_cells[2];
	double ***dss, ***uss, ***vss, ***pss;

	num_cells[0] = mesh->getNumCells()[0];
	num_cells[1] = mesh->getNumCells()[1];

	dss = new double**[2];
	uss = new double**[2];
	vss = new double**[2];
	pss = new double**[2];
	for(size_t i = 0; i < 2; i++)
	{
		dss[i] = new double*[num_cells[0]];
		uss[i] = new double*[num_cells[0]];
		vss[i] = new double*[num_cells[0]];
		pss[i] = new double*[num_cells[0]];
	}
	for(size_t i = 0; i < 2; i++)
		for(size_t j = 0; j < num_cells[0]; j++)
		{
			dss[i][j] = new double[num_cells[1]];
			uss[i][j] = new double[num_cells[1]];
			vss[i][j] = new double[num_cells[1]];
			pss[i][j] = new double[num_cells[1]];
		}


	R = V_init[0];
	U = V_init[1];
	V = V_init[2];
	P = V_init[3];

	boundary_conditions(num_cells, dss, uss, pss, R, U, P);
	linear_solver(num_cells, R, U, V, P, dss, uss, vss, pss);

	double **UFLUX, **VFLUX, **FR, **FRU, **FRV, **FRE;

#pragma omp parallel for collapse(2)
	for(int i = 0; i <= numcells; i++)
	{
		UFLUX[i][j] = uss[0][i][j];
		FR[i] = dss[i] * uss[i];
		FRU[i] = dss[i] * uss[i] * uss[i] + pss[i];
		FRE[i] = (pss[i] / (GAMMA - 1.0) + 0.5*dss[i] * uss[i] * uss[i])*uss[i] + pss[i] * uss[i];
	}
}