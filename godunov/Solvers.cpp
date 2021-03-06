#include "includes\Solvers.h"
#include "includes\flow_parameters.h"

CSolvers::CSolvers()
{
}

CSolvers::~CSolvers(void)
{
}

void CSolvers::getCellsFlowsSecondOrder(int* numcells, CMeshGenerator* mesh, double** R, double** U, double** V, double** P, double*** dss, double*** uss, 
							 double*** vss, double*** pss, boundaryType x_bound_type, boundaryType y_bound_type)
{
	double **C = new double*[numcells[0] + 4];
	double **RC = new double*[numcells[0] + 4];
	double **H = new double*[numcells[0] + 4];

	for(size_t i = 0; i < numcells[0] + 4; i++)
	{
		C[i] = new double[numcells[1] + 4];
		RC[i] = new double[numcells[1] + 4];
		H[i] = new double[numcells[1] + 4];
	}

	CBoundaryConds *boundary = new CBoundaryConds(x_bound_type, y_bound_type);
	double ***flow_data_with_dummy;//R,U,V,P
	double ***R_with_dummy_reconstr, ***U_with_dummy_reconstr, ***V_with_dummy_reconstr, ***P_with_dummy_reconstr;

	R_with_dummy_reconstr = new double**[4];
	U_with_dummy_reconstr = new double**[4];
	V_with_dummy_reconstr = new double**[4];
	P_with_dummy_reconstr = new double**[4];
	for(size_t i = 0; i < 4; i++)
	{
		R_with_dummy_reconstr[i] = new double*[numcells[0] + 4];
		U_with_dummy_reconstr[i] = new double*[numcells[0] + 4];
		V_with_dummy_reconstr[i] = new double*[numcells[0] + 4];
		P_with_dummy_reconstr[i] = new double*[numcells[0] + 4];
	}
	for(size_t i = 0; i < 4; i++)
		for(size_t j = 0; j < numcells[0] + 4; j++)
	{
		R_with_dummy_reconstr[i][j] = new double[numcells[1] + 4];
		U_with_dummy_reconstr[i][j] = new double[numcells[1] + 4];
		V_with_dummy_reconstr[i][j] = new double[numcells[1] + 4];
		P_with_dummy_reconstr[i][j] = new double[numcells[1] + 4];
	}
		

	flow_data_with_dummy = boundary->getFlowDataWithDummy(R, U, V, P, numcells[0], numcells[1]);
	flowDataReconstruction(numcells, flow_data_with_dummy[0], flow_data_with_dummy[1], flow_data_with_dummy[2], 
							flow_data_with_dummy[3], R_with_dummy_reconstr, U_with_dummy_reconstr, 
							V_with_dummy_reconstr, P_with_dummy_reconstr, mesh);

	//for(int k = 0; k < 2; k++)
	//	for(int i = 0; i < numcells[0] + 4; i++)
	//	{
	//		NTracer::traceToFile(U_with_dummy_reconstr[k][i], "double", numcells[1] + 4);
	//	}
	//exit(1);

#pragma omp parallel for collapse(2)
	for(int i = 2; i < numcells[0] + 3; i++)
		for(int j = 2; j < numcells[1] + 3; j++)
		{
			C[i - 1][j] = sqrt(GAMMA*P_with_dummy_reconstr[0][i - 1][j] / R_with_dummy_reconstr[0][i - 1][j]);
			RC[i - 1][j] = R_with_dummy_reconstr[0][i - 1][j] * C[i - 1][j];
			H[i - 1][j] = 1.0 / (RC[i - 1][j]);

			C[i][j] = sqrt(GAMMA*flow_data_with_dummy[3][i][j] / flow_data_with_dummy[0][i][j]);

			C[i][j - 1] = sqrt(GAMMA*P_with_dummy_reconstr[1][i][j - 1] / R_with_dummy_reconstr[1][i][j - 1]);
			RC[i][j - 1] = R_with_dummy_reconstr[1][i][j - 1] * C[i][j - 1];
			H[i][j - 1] = 1.0 / (RC[i][j - 1]);

			RC[i][j] = flow_data_with_dummy[0][i][j] * C[i][j];

			H[i][j] = 1.0 / (RC[i][j]);

			//one direction
			if(U_with_dummy_reconstr[0][i - 1][j] > C[i - 1][j])
			{
				pss[0][i - 2][j - 2] = P_with_dummy_reconstr[0][i - 1][j];
				uss[0][i - 2][j - 2] = U_with_dummy_reconstr[0][i - 1][j];
				dss[0][i - 2][j - 2] = R_with_dummy_reconstr[0][i - 1][j];
				vss[0][i - 2][j - 2] = V_with_dummy_reconstr[0][i - 1][j];
			}
			else if (flow_data_with_dummy[1][i][j] < -C[i][j])
			{
				pss[0][i - 2][j - 2] = flow_data_with_dummy[3][i][j];
				uss[0][i - 2][j - 2] = flow_data_with_dummy[1][i][j];
				dss[0][i - 2][j - 2] = flow_data_with_dummy[0][i][j];
				vss[0][i - 2][j - 2] = flow_data_with_dummy[2][i][j];
			}
			else
			{
				pss[0][i - 2][j - 2] = (U_with_dummy_reconstr[0][i - 1][j] - flow_data_with_dummy[1][i][j] + P_with_dummy_reconstr[0][i - 1][j] * H[i - 1][j] + flow_data_with_dummy[3][i][j] * H[i][j]) / 
					(H[i - 1][j] + H[i][j]);
				uss[0][i - 2][j - 2] = (RC[i - 1][j] * U_with_dummy_reconstr[0][i - 1][j] + RC[i][j] * flow_data_with_dummy[1][i][j] + P_with_dummy_reconstr[0][i - 1][j] - flow_data_with_dummy[3][i][j]) / 
					(RC[i - 1][j] + RC[i][j]);

				if (uss[0][i - 2][j - 2] > 0)
				{
					dss[0][i - 2][j - 2] = R_with_dummy_reconstr[0][i - 1][j] * (1 + (pss[0][i - 2][j - 2] - P_with_dummy_reconstr[0][i - 1][j])/(RC[i - 1][j] * 
						C[i - 1][j]));
					vss[0][i - 2][j - 2] = V_with_dummy_reconstr[0][i - 1][j];
				}
				else
				{
					dss[0][i - 2][j - 2] = flow_data_with_dummy[0][i][j] * (1 + (pss[0][i - 2][j - 2] - flow_data_with_dummy[3][i][j])/(RC[i][j] * C[i][j]));
					vss[0][i - 2][j - 2] = flow_data_with_dummy[2][i][j];
				}
			}

			//second direction
			if(V_with_dummy_reconstr[1][i][j - 1] > C[i][j - 1])
			{
				pss[1][i - 2][j - 2] = P_with_dummy_reconstr[1][i][j - 1];
				uss[1][i - 2][j - 2] = U_with_dummy_reconstr[1][i][j - 1];
				dss[1][i - 2][j - 2] = R_with_dummy_reconstr[1][i][j - 1];
				vss[1][i - 2][j - 2] = V_with_dummy_reconstr[1][i][j - 1];
			}
			else if (flow_data_with_dummy[2][i][j] < -C[i][j])
			{
				pss[1][i - 2][j - 2] = flow_data_with_dummy[3][i][j];
				uss[1][i - 2][j - 2] = flow_data_with_dummy[1][i][j];
				dss[1][i - 2][j - 2] = flow_data_with_dummy[0][i][j];
				vss[1][i - 2][j - 2] = flow_data_with_dummy[2][i][j];
			}
			else
			{
				pss[1][i - 2][j - 2] = (V_with_dummy_reconstr[1][i][j - 1] - flow_data_with_dummy[2][i][j] + P_with_dummy_reconstr[1][i][j - 1] * H[i][j - 1] + flow_data_with_dummy[3][i][j] * H[i][j]) / 
					(H[i][j - 1] + H[i][j]);
				vss[1][i - 2][j - 2] = (RC[i][j - 1] * V_with_dummy_reconstr[1][i][j - 1] + RC[i][j] * flow_data_with_dummy[2][i][j] + P_with_dummy_reconstr[1][i][j - 1] - flow_data_with_dummy[3][i][j]) / 
					(RC[i][j - 1] + RC[i][j]);

				if (vss[1][i - 2][j - 2] > 0)
				{
					dss[1][i - 2][j - 2] = R_with_dummy_reconstr[1][i][j - 1] * (1 + (pss[1][i - 2][j - 2] - P_with_dummy_reconstr[1][i][j - 1])/(RC[i][j - 1] * 
						C[i][j - 1]));
					uss[1][i - 2][j - 2] = U_with_dummy_reconstr[1][i][j - 1];
				}
				else
				{
					dss[1][i - 2][j - 2] = flow_data_with_dummy[0][i][j] * (1 + (pss[1][i - 2][j - 2] - flow_data_with_dummy[3][i][j])/(RC[i][j] * C[i][j]));
					uss[1][i - 2][j - 2] = flow_data_with_dummy[1][i][j];
				}
			}
		}


	for(size_t i = 0; i < 4; i++)
		for(size_t j = 0; j < numcells[0] + 4; j++)
	{
		delete[] R_with_dummy_reconstr[i][j];
		delete[] U_with_dummy_reconstr[i][j];
		delete[] V_with_dummy_reconstr[i][j];
		delete[] P_with_dummy_reconstr[i][j];
	}
	for(size_t i = 0; i < 4; i++)
	{
		delete[] R_with_dummy_reconstr[i];
		delete[] U_with_dummy_reconstr[i];
		delete[] V_with_dummy_reconstr[i];
		delete[] P_with_dummy_reconstr[i];
	}
	delete[] R_with_dummy_reconstr;
	delete[] U_with_dummy_reconstr;
	delete[] V_with_dummy_reconstr;
	delete[] P_with_dummy_reconstr;

	for(size_t i = 0; i < numcells[0] + 4; i++)
	{
		delete[] RC[i];
		delete[] C[i];
		delete[] H[i];
	}

	delete[] RC;
	delete[] C;
	delete[] H;

	delete boundary;

}

void CSolvers::linearSolver(int* numcells, double** R, double** U, double** V, double** P, double*** dss, double*** uss, 
							 double*** vss, double*** pss)
{
	double **C = new double*[numcells[0]];
	double **RC = new double*[numcells[0]];
	double **H = new double*[numcells[0]];

	for(size_t i = 0; i < numcells[0]; i++)
	{
		C[i] = new double[numcells[1]];
		RC[i] = new double[numcells[1]];
		H[i] = new double[numcells[1]];
	}

#pragma omp parallel for collapse(2)
	for(int i = 0; i < numcells[0]; i++)
		for(int j = 0; j < numcells[1]; j++)
		{
			if(i == 0 & j == 0)
				continue;

			if(i > 0)
			{
				C[i - 1][j] = sqrt(GAMMA*P[i - 1][j] / R[i - 1][j]);
				RC[i - 1][j] = R[i - 1][j] * C[i - 1][j];
				H[i - 1][j] = 1.0 / (RC[i - 1][j]);
			}

			C[i][j] = sqrt(GAMMA*P[i][j] / R[i][j]);

			if(j > 0)
			{
				C[i][j - 1] = sqrt(GAMMA*P[i][j - 1] / R[i][j - 1]);
				RC[i][j - 1] = R[i][j - 1] * C[i][j - 1];
				H[i][j - 1] = 1.0 / (RC[i][j - 1]);
			}

			RC[i][j] = R[i][j] * C[i][j];

			H[i][j] = 1.0 / (RC[i][j]);

			//one direction
			if(i > 0)
			{
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
			}

			if(j > 0)
			{
				//second direction
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

					if (vss[1][i][j] > 0)
					{
						dss[1][i][j] = R[i][j - 1] * (1 + (pss[1][i][j] - P[i][j - 1])/(RC[i][j - 1] * 
							C[i][j - 1]));
						uss[1][i][j] = U[i][j - 1];
					}
					else
					{
						dss[1][i][j] = R[i][j] * (1 + (pss[1][i][j] - P[i][j])/(RC[i][j] * C[i][j]));
						uss[1][i][j] = U[i][j];
					}
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

void CSolvers::solve(double*** V_init, int n_steps, int n_save, CMeshGenerator* mesh, boundaryType x_bound_type, boundaryType y_bound_type)
{
	int num_cells[2];
	double ***dss, ***uss, ***vss, ***pss;
	double dx, dy;
	double *y;
	const double G = .01, k = 2;

	CBoundaryConds *bound = new CBoundaryConds(x_bound_type, y_bound_type);
	COutput *output = new COutput;
	CFlowParameters *flow_parameters = new CFlowParameters(mesh);

	num_cells[0] = mesh->getNumCells()[0];
	num_cells[1] = mesh->getNumCells()[1];

	dx = mesh->getMeshStep(0);
	dy = mesh->getMeshStep(1);
	y = mesh->getMeshComponent(1);

	dss = new double**[2];
	uss = new double**[2];
	vss = new double**[2];
	pss = new double**[2];
	for(size_t i = 0; i < 2; i++)
	{
		dss[i] = new double*[num_cells[0] + 1];
		uss[i] = new double*[num_cells[0] + 1];
		vss[i] = new double*[num_cells[0] + 1];
		pss[i] = new double*[num_cells[0] + 1];
	}
	for(size_t i = 0; i < 2; i++)
		for(size_t j = 0; j < num_cells[0] + 1; j++)
		{
			dss[i][j] = new double[num_cells[1] + 1];
			uss[i][j] = new double[num_cells[1] + 1];
			vss[i][j] = new double[num_cells[1] + 1];
			pss[i][j] = new double[num_cells[1] + 1];
		}

	double ***FR, ***FRU, ***FRV, ***FRE;

	FR = new double**[2];
	FRU = new double**[2];
	FRV = new double**[2];
	FRE = new double**[2];
	for(size_t i = 0; i < 2; i++)
	{
		FR[i] = new double*[num_cells[0] + 1];
		FRU[i] = new double*[num_cells[0] + 1];
		FRV[i] = new double*[num_cells[0] + 1];
		FRE[i] = new double*[num_cells[0] + 1];
	}
	for(size_t i = 0; i < 2; i++)
		for(size_t j = 0; j < num_cells[0] + 1; j++)
		{
			FR[i][j] = new double[num_cells[1] + 1];
			FRU[i][j] = new double[num_cells[1] + 1];
			FRV[i][j] = new double[num_cells[1] + 1];
			FRE[i][j] = new double[num_cells[1] + 1];
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

	double **U, **V, **P;

	U = new double*[num_cells[0]];
	V = new double*[num_cells[0]];
	P = new double*[num_cells[0]];
	for(size_t i = 0; i < num_cells[0]; i++)
	{
		U[i] = new double[num_cells[1]];
		V[i] = new double[num_cells[1]];
		P[i] = new double[num_cells[1]];
	}

	R = V_init[0];
	U = V_init[1];
	V = V_init[2];
	P = V_init[3];

#pragma omp parallel for collapse(2)
	for(int i = 0; i < num_cells[0]; i++)
		for(int j = 0; j < num_cells[1]; j++)
	{
		RU[i][j] = R[i][j] * U[i][j];
		RV[i][j] = R[i][j] * V[i][j];
		RE[i][j] = P[i][j] / (GAMMA - 1.0) + 0.5 * R[i][j] * U[i][j] * U[i][j] + 0.5 * R[i][j] * V[i][j] * V[i][j];
	}

	double time = 0.;

	for(int step = 1; step < n_steps; step++)
	{
		double dt = getTimeIncrement(R, U, V, P, mesh);
		time += dt;

		printf("step %d, time %lf\n", step, time);
		//bound->boundary_flows(R, U, V, P, dss, uss, vss, pss, num_cells);
		//linearSolver(num_cells, R, U, V, P, dss, uss, vss, pss);
		getCellsFlowsSecondOrder(num_cells, mesh, R, U, V, P, dss, uss, vss, pss, x_bound_type, y_bound_type);

		//double *traceData = new double[num_cells[1]];
		//for(int k = 0; k < 2; k++)
		//	for(int i = 0; i < num_cells[0] + 1; i++)
		//	{
		//		NTracer::traceToFile(pss[k][i], "double", num_cells[1]);
		//	}
		//exit(1);

#pragma omp parallel for collapse(2)
		for(int i = 0; i <= num_cells[0]; i++)
			for(int j = 0; j <= num_cells[1]; j++)
		{
			if(j < num_cells[1])
			{
				FR[0][i][j] = dss[0][i][j] * uss[0][i][j];
				FRU[0][i][j] = dss[0][i][j] * uss[0][i][j] * uss[0][i][j] + pss[0][i][j];
				FRV[0][i][j] = dss[0][i][j] * uss[0][i][j] * vss[0][i][j];
				FRE[0][i][j] = (pss[0][i][j] / (GAMMA - 1.0) + 0.5*dss[0][i][j] * (uss[0][i][j] * uss[0][i][j] + vss[0][i][j] * vss[0][i][j]))*uss[0][i][j] + 
					pss[0][i][j] * uss[0][i][j];
			}

			if(i < num_cells[0])
			{
				FR[1][i][j] = dss[1][i][j] * vss[1][i][j];
				FRV[1][i][j] = dss[1][i][j] * vss[1][i][j] * vss[1][i][j] + pss[1][i][j];
				FRU[1][i][j] = dss[1][i][j] * uss[1][i][j] * vss[1][i][j];
				FRE[1][i][j] = (pss[1][i][j] / (GAMMA - 1.0) + 0.5*dss[1][i][j] * (uss[1][i][j] * uss[1][i][j] + vss[1][i][j] * vss[1][i][j]))*vss[1][i][j] + 
					pss[1][i][j] * vss[1][i][j];
			}
		}
		//delete []traceData;
		//exit(1);

#pragma omp parallel for collapse(2)
		for(int i = 0; i < num_cells[0]; i++)
			for(int j = 0; j < num_cells[1]; j++)
		{
			R[i][j] = R[i][j] - dt/dx * (FR[0][i + 1][j] - FR[0][i][j]) - dt/dy * (FR[1][i][j + 1] - FR[1][i][j]);
			RU[i][j] = RU[i][j] - dt/dx * (FRU[0][i + 1][j] - FRU[0][i][j]) - dt/dy * (FRU[1][i][j + 1] - FRU[1][i][j]) + dt * G * R[i][j] * sin(k * y[j]);
			RV[i][j] = RV[i][j] - dt/dx * (FRV[0][i + 1][j] - FRV[0][i][j]) - dt/dy * (FRV[1][i][j + 1] - FRV[1][i][j]);// - dt * R[i][j] * .1;
			RE[i][j] = RE[i][j] - dt/dx * (FRE[0][i + 1][j] - FRE[0][i][j]) - dt/dy * (FRE[1][i][j + 1] - FRE[1][i][j]);
		}

#pragma omp parallel for collapse(2)
		for(int i = 0; i < num_cells[0]; i++)
			for(int j = 0; j < num_cells[1]; j++)
			{
				U[i][j] = RU[i][j] / R[i][j];
				V[i][j] = RV[i][j] / R[i][j];
				P[i][j] = (GAMMA - 1.0) * (RE[i][j] - 0.5 * RU[i][j] * U[i][j] - 0.5 * RV[i][j] * V[i][j]);

				if(P[i][j] < 0)
				{
					std::cout << "negative pressure at " << "i " << i << ", j " << j;
					exit(1);
				}

				if(R[i][j] < 0)
				{
					std::cout << "negative density at " << "i " << i << ", j " << j;
					exit(1);
				}
			}

		if((step % n_save) == 0/*time > .1*/)
		{
			std::stringstream namefile_2d /*namefile_x, namefile_y*/;
			//namefile_x << "x_" << step << ".dat";
			//namefile_y << "y_" << step << ".dat";
			namefile_2d << "2d_" << step << ".dat";
			//output->save1dPlotXAxis(namefile_x.str().c_str(), mesh, time, 50, R, U, V, P);
			//output->save1dPlotYAxis(namefile_y.str().c_str(), mesh, time, 50, R, U, V, P);
			double **S = flow_parameters->getFlowEntropy(R, P);
			double **S_diff = flow_parameters->getFlowEntropyDiff();
			output->save2dPlot(namefile_2d.str().c_str(), mesh, time, R, U, V, P, S, S_diff);
			//exit(1);
		}
		//flow_parameters->saveFlowEntropyCurrTimeStep(R, P);

		//if((step % 10) == 0)
		//{
		//	output->saveKineticEnergyAndEnstrophy("kolm.dat", mesh, R, U, V, time);
		//}

		//const int n_save_flow_parameters = 10000;

		//if((step % n_save_flow_parameters) == 0 || step == 1)
		//{
			//std::stringstream namefile_spectrum;
			//double* E_k;

			//namefile_spectrum << "flow_spectrum_" << time << ".dat";
			//E_k = flow_parameters->getKineticEnergySpectrum(U, V);
			//for(int k = 0; k < flow_parameters->k_max(); k++)
			//{
			//	output->saveKineticEnergySpectrum(namefile_spectrum.str().c_str(), E_k[k], k + 1);
			//}
		//}
		output->saveData("flow_energy_power.dat", flow_parameters->getEnergyPower(U, R, G, k), time);
		output->saveData("flow_enstrophy_power.dat", flow_parameters->getEnstrophyPower(U, V, R, G, k), time);
	}

	delete bound;
	delete output;
	delete flow_parameters;

	for(size_t i = 0; i < 2; i++)
		for(size_t j = 0; j < num_cells[0] + 1; j++)
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
		for(size_t j = 0; j < num_cells[0] + 1; j++)
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

	for(size_t i = 0; i < num_cells[0]; i++)
	{
		delete []R[i];
		delete []RU[i];
		delete []RV[i];
		delete []RE[i];
	}
	delete []R;
	delete []RU;
	delete []RV;
	delete []RE;


	for(size_t i = 0; i < num_cells[0]; i++)
	{
		delete []U[i];
		delete []V[i];
		delete []P[i];
	}
	delete []U;
	delete []V;
	delete []P;
}

double CSolvers::getTimeIncrement(double** R, double** U, double** V, double** P, CMeshGenerator* mesh)
{
	double c;
	double velocity_max = U[0][0];

    for(int i = 0; i < mesh->getNumCells()[0]; i++)
    {
		for(int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			c = sqrt(GAMMA * P[i][j]/R[i][j]);
			velocity_max = U[i][j] > velocity_max ? U[i][j] : velocity_max;
			velocity_max = V[i][j] > velocity_max ? V[i][j] : velocity_max;
			velocity_max = V[i][j] - c > velocity_max ? V[i][j] - c : velocity_max;
			velocity_max = V[i][j] + c > velocity_max ? V[i][j] + c : velocity_max;
			velocity_max = U[i][j] - c > velocity_max ? U[i][j] - c : velocity_max;
			velocity_max = U[i][j] + c > velocity_max ? U[i][j] + c : velocity_max;

		}
	}

	return CFL * std::min(mesh->getMeshStep(0), mesh->getMeshStep(1)) / velocity_max;
}

void CSolvers::flowDataReconstruction(int* numcells, double** R_with_dummy, double** U_with_dummy, double** V_with_dummy, 
							double** P_with_dummy, double*** R_with_dummy_reconstr, double*** U_with_dummy_reconstr, 
							double*** V_with_dummy_reconstr, double*** P_with_dummy_reconstr, CMeshGenerator* mesh)
{
	double ***R_with_dummy_grad, ***U_with_dummy_grad, ***V_with_dummy_grad, ***P_with_dummy_grad;
	double ***R_limiter, ***U_limiter, ***V_limiter, ***P_limiter;
	int x_cells_num = mesh->getNumCells()[0];
	int y_cells_num = mesh->getNumCells()[1];
	double dx = mesh->getMeshStep(0);
	double dy = mesh->getMeshStep(1);

	R_with_dummy_grad = new double**[2];
	U_with_dummy_grad = new double**[2];
	V_with_dummy_grad = new double**[2];
	P_with_dummy_grad = new double**[2];
	for(int i = 0; i < 2; i++)
	{
		R_with_dummy_grad[i] = new double*[x_cells_num + 4];
		U_with_dummy_grad[i] = new double*[x_cells_num + 4];
		V_with_dummy_grad[i] = new double*[x_cells_num + 4];
		P_with_dummy_grad[i] = new double*[x_cells_num + 4];
	}
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < x_cells_num + 4; j++)
		{
			R_with_dummy_grad[i][j] = new double[y_cells_num + 4];
			U_with_dummy_grad[i][j] = new double[y_cells_num + 4];
			V_with_dummy_grad[i][j] = new double[y_cells_num + 4];
			P_with_dummy_grad[i][j] = new double[y_cells_num + 4];
		}
	}

	R_limiter = new double**[2];
	U_limiter = new double**[2];
	V_limiter = new double**[2];
	P_limiter = new double**[2];
	for(int i = 0; i < 2; i++)
	{
		R_limiter[i] = new double*[x_cells_num + 2];
		U_limiter[i] = new double*[x_cells_num + 2];
		V_limiter[i] = new double*[x_cells_num + 2];
		P_limiter[i] = new double*[x_cells_num + 2];
	}
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < x_cells_num + 2; j++)
		{
			R_limiter[i][j] = new double[y_cells_num + 2];
			U_limiter[i][j] = new double[y_cells_num + 2];
			V_limiter[i][j] = new double[y_cells_num + 2];
			P_limiter[i][j] = new double[y_cells_num + 2];
		}
	}
	getFlowDataGradient(R_with_dummy, U_with_dummy, V_with_dummy, P_with_dummy, R_with_dummy_grad, 
						U_with_dummy_grad, V_with_dummy_grad, P_with_dummy_grad, mesh);
	flowVariablesLimiters(R_with_dummy, U_with_dummy, V_with_dummy, P_with_dummy, R_limiter, U_limiter, 
						V_limiter, P_limiter, mesh);

#pragma omp parallel for collapse(2)
	for(int i = 1; i < x_cells_num + 3; i++)
	{
		for(int j = 1; j < y_cells_num + 3; j++)
		{
			U_with_dummy_reconstr[0][i][j] = U_with_dummy[i][j] - .5 * dx * U_limiter[0][i - 1][j - 1] * U_with_dummy_grad[0][i][j];
			U_with_dummy_reconstr[1][i][j] = U_with_dummy[i][j] + .5 * dy * U_limiter[1][i - 1][j - 1] * U_with_dummy_grad[1][i][j];
			U_with_dummy_reconstr[2][i][j] = U_with_dummy[i][j] + .5 * dx * U_limiter[0][i - 1][j - 1] * U_with_dummy_grad[0][i][j];
			U_with_dummy_reconstr[3][i][j] = U_with_dummy[i][j] - .5 * dy * U_limiter[1][i - 1][j - 1] * U_with_dummy_grad[1][i][j];

			V_with_dummy_reconstr[0][i][j] = V_with_dummy[i][j] - .5 * dx * V_limiter[0][i - 1][j - 1] * V_with_dummy_grad[0][i][j];
			V_with_dummy_reconstr[1][i][j] = V_with_dummy[i][j] + .5 * dy * V_limiter[1][i - 1][j - 1] * V_with_dummy_grad[1][i][j];
			V_with_dummy_reconstr[2][i][j] = V_with_dummy[i][j] + .5 * dx * V_limiter[0][i - 1][j - 1] * V_with_dummy_grad[0][i][j];
			V_with_dummy_reconstr[3][i][j] = V_with_dummy[i][j] - .5 * dy * V_limiter[1][i - 1][j - 1] * V_with_dummy_grad[1][i][j];

			R_with_dummy_reconstr[0][i][j] = R_with_dummy[i][j] - .5 * dx * R_limiter[0][i - 1][j - 1] * R_with_dummy_grad[0][i][j];
			R_with_dummy_reconstr[1][i][j] = R_with_dummy[i][j] + .5 * dy * R_limiter[1][i - 1][j - 1] * R_with_dummy_grad[1][i][j];
			R_with_dummy_reconstr[2][i][j] = R_with_dummy[i][j] + .5 * dx * R_limiter[0][i - 1][j - 1] * R_with_dummy_grad[0][i][j];
			R_with_dummy_reconstr[3][i][j] = R_with_dummy[i][j] - .5 * dy * R_limiter[1][i - 1][j - 1] * R_with_dummy_grad[1][i][j];

			P_with_dummy_reconstr[0][i][j] = P_with_dummy[i][j] - .5 * dx * P_limiter[0][i - 1][j - 1] * P_with_dummy_grad[0][i][j];
			P_with_dummy_reconstr[1][i][j] = P_with_dummy[i][j] + .5 * dy * P_limiter[1][i - 1][j - 1] * P_with_dummy_grad[1][i][j];
			P_with_dummy_reconstr[2][i][j] = P_with_dummy[i][j] + .5 * dx * P_limiter[0][i - 1][j - 1] * P_with_dummy_grad[0][i][j];
			P_with_dummy_reconstr[3][i][j] = P_with_dummy[i][j] - .5 * dy * P_limiter[1][i - 1][j - 1] * P_with_dummy_grad[1][i][j];
		}
	}


	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < x_cells_num + 4; j++)
		{
			delete[] R_with_dummy_grad[i][j];
			delete[] U_with_dummy_grad[i][j];
			delete[] V_with_dummy_grad[i][j];
			delete[] P_with_dummy_grad[i][j];
		}
	}
	for(int i = 0; i < 2; i++)
	{
		delete[] R_with_dummy_grad[i];
		delete[] U_with_dummy_grad[i];
		delete[] V_with_dummy_grad[i];
		delete[] P_with_dummy_grad[i];
	}
	delete[] R_with_dummy_grad;
	delete[] U_with_dummy_grad;
	delete[] V_with_dummy_grad;
	delete[] P_with_dummy_grad;

	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < x_cells_num + 2; j++)
		{
			delete[] R_limiter[i][j];
			delete[] U_limiter[i][j];
			delete[] V_limiter[i][j];
			delete[] P_limiter[i][j];
		}
	}
	for(int i = 0; i < 2; i++)
	{
		delete[] R_limiter[i];
		delete[] U_limiter[i];
		delete[] V_limiter[i];
		delete[] P_limiter[i];
	}
	delete[] R_limiter;
	delete[] U_limiter;
	delete[] V_limiter;
	delete[] P_limiter;
}

void CSolvers::getFlowDataGradient(double** R_with_dummy, double** U_with_dummy, double** V_with_dummy, 
						double** P_with_dummy, double*** R_with_dummy_grad, double*** U_with_dummy_grad, 
						double*** V_with_dummy_grad, double*** P_with_dummy_grad, CMeshGenerator* mesh)
{
	int x_cells_num = mesh->getNumCells()[0];
	int y_cells_num = mesh->getNumCells()[1];
	double dx = mesh->getMeshStep(0);
	double dy = mesh->getMeshStep(1);

#pragma omp parallel for collapse(2)
	for(int i = 1; i < x_cells_num + 3; i++)
	{
		for(int j = 1; j < y_cells_num + 3; j++)
		{
			R_with_dummy_grad[0][i][j] = .5 * (R_with_dummy[i + 1][j] - R_with_dummy[i - 1][j])/dx;
			R_with_dummy_grad[1][i][j] = .5 * (R_with_dummy[i][j + 1] - R_with_dummy[i][j - 1])/dy;
			U_with_dummy_grad[0][i][j] = .5 * (U_with_dummy[i + 1][j] - U_with_dummy[i - 1][j])/dx;
			U_with_dummy_grad[1][i][j] = .5 * (U_with_dummy[i][j + 1] - U_with_dummy[i][j - 1])/dy;
			V_with_dummy_grad[0][i][j] = .5 * (V_with_dummy[i + 1][j] - V_with_dummy[i - 1][j])/dx;
			V_with_dummy_grad[1][i][j] = .5 * (V_with_dummy[i][j + 1] - V_with_dummy[i][j - 1])/dy;
			P_with_dummy_grad[0][i][j] = .5 * (P_with_dummy[i + 1][j] - P_with_dummy[i - 1][j])/dx;
			P_with_dummy_grad[1][i][j] = .5 * (P_with_dummy[i][j + 1] - P_with_dummy[i][j - 1])/dy;
		}
	}
}

void CSolvers::flowVariablesLimiters(double** R_with_dummy, double** U_with_dummy, double** V_with_dummy, double** P_with_dummy, 
											  double*** R_limiter, double*** U_limiter, double*** V_limiter, double*** P_limiter, 
											  CMeshGenerator* mesh)
{
	int x_cells_num = mesh->getNumCells()[0];
	int y_cells_num = mesh->getNumCells()[1];

#pragma omp parallel for collapse(2)
	for(int i = 1; i < x_cells_num + 3; i++)
	{
		for(int j = 1; j < y_cells_num + 3; j++)
		{
			R_limiter[0][i - 1][j - 1] = std::max((double)0, std::min((double)1, (R_with_dummy[i][j] - R_with_dummy[i - 1][j])/(R_with_dummy[i + 1][j] - R_with_dummy[i][j])));
			R_limiter[1][i - 1][j - 1] = std::max((double)0, std::min((double)1, (R_with_dummy[i][j] - R_with_dummy[i][j - 1])/(R_with_dummy[i][j + 1] - R_with_dummy[i][j])));
			U_limiter[0][i - 1][j - 1] = std::max((double)0, std::min((double)1, (U_with_dummy[i][j] - U_with_dummy[i - 1][j])/(U_with_dummy[i + 1][j] - U_with_dummy[i][j])));
			U_limiter[1][i - 1][j - 1] = std::max((double)0, std::min((double)1, (U_with_dummy[i][j] - U_with_dummy[i][j - 1])/(U_with_dummy[i][j + 1] - U_with_dummy[i][j])));
			V_limiter[0][i - 1][j - 1] = std::max((double)0, std::min((double)1, (V_with_dummy[i][j] - V_with_dummy[i - 1][j])/(V_with_dummy[i + 1][j] - V_with_dummy[i][j])));
			V_limiter[1][i - 1][j - 1] = std::max((double)0, std::min((double)1, (V_with_dummy[i][j] - V_with_dummy[i][j - 1])/(V_with_dummy[i][j + 1] - V_with_dummy[i][j])));
			P_limiter[0][i - 1][j - 1] = std::max((double)0, std::min((double)1, (P_with_dummy[i][j] - P_with_dummy[i - 1][j])/(P_with_dummy[i + 1][j] - P_with_dummy[i][j])));
			P_limiter[1][i - 1][j - 1] = std::max((double)0, std::min((double)1, (P_with_dummy[i][j] - P_with_dummy[i][j - 1])/(P_with_dummy[i][j + 1] - P_with_dummy[i][j])));
		}
	}
}