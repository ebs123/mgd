#include "includes\Solvers.h"
#include "includes\flow_parameters.h"

CSolvers::CSolvers()
{
}

CSolvers::~CSolvers(void)
{
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

void CSolvers::solve(double*** V_init, int n_steps, int n_save, CMeshGenerator* mesh)
{
	int num_cells[2];
	double ***dss, ***uss, ***vss, ***pss;
	double dx, dy;
	double *y;
	const double G = 0.01, k = 1;

	CBoundaryConds *bound = new CBoundaryConds;
	COutput *output = new COutput;

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

	double **U, **V, **P, **S;

	U = new double*[num_cells[0]];
	V = new double*[num_cells[0]];
	P = new double*[num_cells[0]];
	S = new double*[num_cells[0]];
	for(size_t i = 0; i < num_cells[0]; i++)
	{
		U[i] = new double[num_cells[1]];
		V[i] = new double[num_cells[1]];
		P[i] = new double[num_cells[1]];
		S[i] = new double[num_cells[1]];
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
		bound->boundary_flows(R, U, V, P, dss, uss, vss, pss, num_cells);
		linearSolver(num_cells, R, U, V, P, dss, uss, vss, pss);

			//double *traceData = new double[10];
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
			RV[i][j] = RV[i][j] - dt/dx * (FRV[0][i + 1][j] - FRV[0][i][j]) - dt/dy * (FRV[1][i][j + 1] - FRV[1][i][j]);
			RE[i][j] = RE[i][j] - dt/dx * (FRE[0][i + 1][j] - FRE[0][i][j]) - dt/dy * (FRE[1][i][j + 1] - FRE[1][i][j]);
		}


#pragma omp parallel for collapse(2)
		for(int i = 0; i < num_cells[0]; i++)
			for(int j = 0; j < num_cells[1]; j++)
			{
				U[i][j] = RU[i][j] / R[i][j];
				V[i][j] = RV[i][j] / R[i][j];
				P[i][j] = (GAMMA - 1.0) * (RE[i][j] - 0.5 * RU[i][j] * U[i][j] - 0.5 * RV[i][j] * V[i][j]);
				S[i][j] = log(P[i][j] / pow(R[i][j], GAMMA));

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

			//traceData[0] = i;
			//traceData[1] = j;
			//traceData[2] = dss[0][i][j];
			//traceData[3] = dss[1][i][j];
			//traceData[4] = uss[0][i][j];
			//traceData[5] = uss[1][i][j];
			//traceData[6] = vss[0][i][j];
			//traceData[7] = vss[1][i][j];
			//traceData[8] = pss[0][i][j];
			//traceData[9] = pss[1][i][j];
			//NTracer::traceToFile(traceData, "double", 10);

			}

		if((step % n_save) == 0)
		{
			std::stringstream namefile_2d;
			/*namefile_x << "x_" << step << ".dat";
			namefile_y << "y_" << step << ".dat";*/
			namefile_2d << "2d_" << step << ".dat";
			/*output->save1dPlotXAxis(namefile_x.str().c_str(), mesh, time, 200, R, U, V, P, P);
			output->save1dPlotYAxis(namefile_y.str().c_str(), mesh, time, 200, R, U, V, P, P);*/
			output->save2dPlot(namefile_2d.str().c_str(), mesh, time, R, U, V, P, P);
		}

		int n_save_flow_parameters = 10000;
		if((step % 10) == 0)
		{
			output->saveKineticEnergyAndEnstrophy("kolm.dat", mesh, R, U, V, time);
		}

		if((step % n_save_flow_parameters) == 0)
		{
			std::stringstream namefile_spectrum;
			CFlowParameters *flow_parameters = new CFlowParameters;
			double E_spectrum;

			namefile_spectrum << "kolm_spectrum_" << time << ".dat";
			for(int k_x = 0; k_x < 110; k_x++)
			{
				int k_y = 1;
				E_spectrum = flow_parameters->getKineticEnergySpectrum((double)k_x, (double)k_y, U, V, mesh);
				output->saveKineticEnergySpectrum(namefile_spectrum.str().c_str(), E_spectrum, (double)k_x, (double)k_y);
			}

			delete flow_parameters;
		}
	}

	delete bound;
	delete output;

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
		delete []S[i];
	}
	delete []U;
	delete []V;
	delete []P;
	delete []S;
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