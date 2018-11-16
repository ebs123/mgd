#include "includes\Solvers.h"

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

void CSolvers::solve(double*** V_init, int n_steps, int n_save, CMeshGenerator* mesh)
{
	int num_cells[2];
	double ***dss, ***uss, ***vss, ***pss;
	double dx, dy;

	CBoundaryConds *bound = new CBoundaryConds;
	COutput *output = new COutput;

	num_cells[0] = mesh->getNumCells()[0];
	num_cells[1] = mesh->getNumCells()[1];

	dx = mesh->getMeshStep(0);
	dy = mesh->getMeshStep(1);

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

	double ***FR, ***FRU, ***FRV, ***FRE;

	FR = new double**[2];
	FRU = new double**[2];
	FRV = new double**[2];
	FRE = new double**[2];
	for(size_t i = 0; i < 2; i++)
	{
		FR[i] = new double*[num_cells[0]];
		FRU[i] = new double*[num_cells[0]];
		FRV[i] = new double*[num_cells[0]];
		FRE[i] = new double*[num_cells[0]];
	}
	for(size_t i = 0; i < 2; i++)
		for(size_t j = 0; j < num_cells[0]; j++)
		{
			FR[i][j] = new double[num_cells[1]];
			FRU[i][j] = new double[num_cells[1]];
			FRV[i][j] = new double[num_cells[1]];
			FRE[i][j] = new double[num_cells[1]];
		}

	double **R, **RU, **RV, **RE;

	R = new double*[num_cells[0]];
	RU = new double*[num_cells[0]];
	RV = new double*[num_cells[0]];
	RE = new double*[num_cells[0]];
	for(size_t i = 0; i < num_cells[0]; i++)
	{
		R[i] = new double[num_cells[1]];
		RU[i] = new double[num_cells[1]];
		RV[i] = new double[num_cells[1]];
		RE[i] = new double[num_cells[1]];
	}

	double **R, **U, **V, **P, **S;

	R = new double*[num_cells[0]];
	U = new double*[num_cells[0]];
	V = new double*[num_cells[0]];
	P = new double*[num_cells[0]];
	S = new double*[num_cells[0]];
	for(size_t i = 0; i < num_cells[0]; i++)
	{
		R[i] = new double[num_cells[1]];
		U[i] = new double[num_cells[1]];
		V[i] = new double[num_cells[1]];
		P[i] = new double[num_cells[1]];
		S[i] = new double[num_cells[1]];
	}

	R = V_init[0];
	U = V_init[1];
	V = V_init[2];
	P = V_init[3];

	for(int step = 1; step < n_steps; step++)
	{
		bound->boundary_flows(R, U, V, P, dss, uss, vss, pss, numcells);
		linear_solver(num_cells, R, U, V, P, dss, uss, vss, pss);

#pragma omp parallel for collapse(2)
			for(int i = 0; i <= numcells[0]; i++)
				for(int j = 0; j <= numcells[1]; j++)
			{
				FR[0][i][j] = dss[0][i][j] * uss[0][i][j];
				FRU[0][i][j] = dss[0][i][j] * uss[0][i][j] * uss[0][i][j] + pss[0][i][j];
				FRV[0][i][j] = dss[0][i][j] * uss[0][i][j] * vss[0][i][j];
				FRE[0][i][j] = (pss[0][i][j] / (GAMMA - 1.0) + 0.5*dss[0][i][j] * (uss[0][i][j] * uss[0][i][j] + vss[0][i][j] * vss[0][i][j]))*uss[0][i][j] + 
					pss[0][i][j] * uss[0][i][j];

				FR[1][i][j] = dss[1][i][j] * vss[1][i][j];
				FRV[1][i][j] = dss[1][i][j] * vss[1][i][j] * vss[1][i][j] + pss[1][i][j];
				FRU[1][i][j] = dss[1][i][j] * uss[1][i][j] * vss[1][i][j];
				FRE[1][i][j] = (pss[1][i][j] / (GAMMA - 1.0) + 0.5*dss[1][i][j] * (uss[1][i][j] * uss[1][i][j] + vss[1][i][j] * vss[1][i][j]))*vss[1][i][j] + 
					pss[1][i][j] * vss[1][i][j];
			}

			double dt = getTimeIncrement(U, V, mesh);

#pragma omp parallel for collapse(2)
			for(int i = 0; i <= numcells[0]; i++)
				for(int j = 0; j <= numcells[1]; j++)
			{
				R[i][j] = R[i][j] - dt/dx * (FR[0][i + 1][j] - FR[0][i][j]) - dt/dy * (FR[1][i + 1][j] - FR[1][i][j]);
				RU[i][j] = RU[i][j] - dt/dx * (FRU[0][i + 1][j] - FRU[0][i][j]) - dt/dy * (FRU[1][i + 1][j] - FRU[1][i][j]);
				RV[i][j] = RV[i][j] - dt/dx * (FRU[0][i + 1][j] - FRU[0][i][j]) - dt/dy * (FRU[1][i + 1][j] - FRU[1][i][j]);
				RE[i][j] = RE[i][j] - dt/dx * (FRE[0][i + 1][j] - FRE[0][i][j]) - dt/dy * (FRE[1][i + 1][j] - FRE[1][i][j]);
			}

#pragma omp parallel for collapse(2)
			for(int i = 0; i <= numcells[0]; i++)
				for(int j = 0; j <= numcells[1]; j++)
				{
					U[i][j] = RU[i][j] / R[i][j];
					V[i][j] = RV[i][j] / R[i][j];
					P[i][j] = (GAMMA - 1.0) * (RE[i][j] - 0.5 * RU[i][j] * U[i][j] - 0.5 * RV[i][j] * V[i][j]);
					S[i][j] = log(P[i][j] / pow(R[i][j], GAMMA));
				}

			if((step % n_save) == 0)
			{
				std::stringstream namefile;
				namefile << step << ".dat";
				output->save2dPlot(namefile.str().c_str(), mesh, R, U, V, P, S);
			}
	}

	delete bound;
	delete output;

	for(size_t i = 0; i < 2; i++)
		for(size_t j = 0; j < num_cells[0]; j++)
		{
			delete []dss[i][j];
			delete []uss[i][j];
			delete []vss[i][j];
			delete []pss[i][j];
		}
	for(size_t i = 0; i < 2; i++)
	{
		delete []dss[i];
		delete []uss[i];
		delete []vss[i];
		delete []pss[i];
	}
	delete []dss;
	delete []uss;
	delete []vss;
	delete []pss;

	for(size_t i = 0; i < 2; i++)
		for(size_t j = 0; j < num_cells[0]; j++)
		{
			delete []FR[i][j];
			delete []FRU[i][j];
			delete []FRV[i][j];
			delete []FRE[i][j];
		}
	for(size_t i = 0; i < 2; i++)
	{
		delete []FR[i];
		delete []FRU[i];
		delete []FRV[i];
		delete []FRE[i];
	}
	delete []FR;
	delete []FRU;
	delete []FRV;
	delete []FRE;
}

double CSolvers::getTimeIncrement(double** U, double** V, CMeshGenerator* mesh)
{
	double velocity_max = U[0][0];

    for(int i = 0; i < mesh->getNumCells()[0]; i++)
    {
		for(int j = 0; j < mesh->getNumCells()[1]; j++)
		{
		  if(U[i][j] > velocity_max)
		  {
			  velocity_max = U[i][j];
		  }

		  if(V[i][j] > velocity_max)
		  {
			  velocity_max = V[i][j];
		  }
		}
	}

	return CFL * min(mesh->getMeshStep(0), mesh->getMeshStep(1)) / velocity_max;
}