#include "includes\flow_parameters.h"

CFlowParameters::CFlowParameters(CMeshGenerator *grid)
{
	k_x_max = 40;
	k_y_max = 40;
	E = new double[(int)sqrt((double)(k_x_max * k_x_max + k_y_max * k_y_max))];
	S = new double*[grid->getNumCells()[0]];
	for(size_t i = 0; i < grid->getNumCells()[0]; i++)
	{
		S[i] = new double[grid->getNumCells()[1]];
	}

	S_last_step = new double*[grid->getNumCells()[0]];
	for(size_t i = 0; i < grid->getNumCells()[0]; i++)
	{
		S_last_step[i] = new double[grid->getNumCells()[1]];
	}

	S_diff = new double*[grid->getNumCells()[0]];
	for(size_t i = 0; i < grid->getNumCells()[0]; i++)
	{
		S_diff[i] = new double[grid->getNumCells()[1]];
	}
	mesh = grid;
}

CFlowParameters::~CFlowParameters(void)
{
	delete E;
	for(size_t i = 0; i < mesh->getNumCells()[0]; i++)
	{
		delete []S[i];
	}
	delete []S;

	for(size_t i = 0; i < mesh->getNumCells()[0]; i++)
	{
		delete []S_diff[i];
	}
	delete []S_diff;

	for(size_t i = 0; i < mesh->getNumCells()[0]; i++)
	{
		delete []S_last_step[i];
	}
	delete []S_last_step;
}

double** CFlowParameters::getFlowEntropy(double** R, double** P)
{
	for(size_t i = 0; i < mesh->getNumCells()[0]; i++)
	{
		for(size_t j = 0; j < mesh->getNumCells()[1]; j++)
		{
			S[i][j] = log(P[i][j] / pow(R[i][j], GAMMA));
		}
	}

	return S;
}

double** CFlowParameters::getFlowEntropyDiff()
{
	for(size_t i = 0; i < mesh->getNumCells()[0]; i++)
	{
		for(size_t j = 0; j < mesh->getNumCells()[1]; j++)
		{
			S_diff[i][j] = S[i][j] - S_last_step[i][j];
		}
	}

	return S_diff;
}

void CFlowParameters::saveFlowEntropyCurrTimeStep(double** R, double** P)
{
	getFlowEntropy(R, P);
	for(size_t i = 0; i < mesh->getNumCells()[0]; i++)
	{
		for(size_t j = 0; j < mesh->getNumCells()[1]; j++)
		{
			S_last_step[i][j] = S[i][j];
		}
	}
}

double* CFlowParameters::getKineticEnergySpectrum(double** U, double** V)
{
	double U_fourier_i[4], V_fourier_i[4];
	const int N = mesh->getNumCells()[0] - 2;
	double h, E_kx_ky;
	const int k_max = sqrt((double)(k_x_max * k_x_max + k_y_max * k_y_max));

	for(int i = 0; i < k_max; i++)
	{
		E[i] = 0;
	}

	h = 2 * pi/(N + 1);

	for(int k_loop = 0; k_loop < k_max; k_loop++)
	{
#pragma omp parallel for collapse(2)
		for(int k_x = 1; k_x < k_x_max + 1; k_x++)
		{
			for(int k_y = 1; k_y < k_y_max + 1; k_y++)
			{
				double k_mod = sqrt((double)(k_x * k_x + k_y * k_y));

				if(abs(k_loop - k_mod + 1) < .5)
				{
					for(int i = 0; i < 4; i++)
					{
						U_fourier_i[i] = 0;
						V_fourier_i[i] = 0;
					}

					#pragma omp parallel for collapse(2)
						for(int i = 1; i < N + 2; i++)
						{
							for(int j = 1; j < N + 2; j++)
							{
								U_fourier_i[0] += 1/pi * U[i][j] * cos(k_x * h * (i - .5)) * cos(k_y * h * (j - .5)) * h * h;
								U_fourier_i[1] += 1/pi * U[i][j] * cos(k_x * h * (i - .5)) * sin(k_y * h * (j - .5)) * h * h;
								U_fourier_i[2] += 1/pi * U[i][j] * sin(k_x * h * (i - .5)) * cos(k_y * h * (j - .5)) * h * h;
								U_fourier_i[3] += 1/pi * U[i][j] * sin(k_x * h * (i - .5)) * sin(k_y * h * (j - .5)) * h * h;

								V_fourier_i[0] += 1/pi * V[i][j] * cos(k_x * h * (i - .5)) * cos(k_y * h * (j - .5)) * h * h;
								V_fourier_i[1] += 1/pi * V[i][j] * cos(k_x * h * (i - .5)) * sin(k_y * h * (j - .5)) * h * h;
								V_fourier_i[2] += 1/pi * V[i][j] * sin(k_x * h * (i - .5)) * cos(k_y * h * (j - .5)) * h * h;
								V_fourier_i[3] += 1/pi * V[i][j] * sin(k_x * h * (i - .5)) * sin(k_y * h * (j - .5)) * h * h;
							}
						}

						E_kx_ky = (U_fourier_i[0] * U_fourier_i[0] + U_fourier_i[1] * U_fourier_i[1] + U_fourier_i[2] * U_fourier_i[2] + 
							U_fourier_i[3] * U_fourier_i[3] + V_fourier_i[0] * V_fourier_i[0] + V_fourier_i[1] * V_fourier_i[1] + 
							V_fourier_i[2] * V_fourier_i[2] + V_fourier_i[3] * V_fourier_i[3]) * .5;

					E[k_loop] += E_kx_ky;
				}
			}
		}
	}

	return E;
}

double CFlowParameters::getEnergyPower(double** U, double** R, const double G, const double k)
{
	double energy_power = 0;
	const int N = mesh->getNumCells()[0] - 2;
	double h;
	double *y;

	y = mesh->getMeshComponent(1);
	h = 2 * pi/(N + 1);

	for (int i = 0; i < mesh->getNumCells()[0]; i++)
	{
		for (int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			energy_power += h * h * G * R[i][j] * sin(k * y[j]) * U[i][j];
		}
	}

	return energy_power;
}

double CFlowParameters::getEnstrophyPower(double** U, double** V, double** R, const double G, const double k)
{
	double enstrophy_power = 0;
	const int N = mesh->getNumCells()[0] - 2;
	double h;
	double *y;
	int j_plus, i_plus;
	int j_minus, i_minus;

	y = mesh->getMeshComponent(1);
	h = 2 * pi/(N + 1);

	for (int i = 0; i < mesh->getNumCells()[0]; i++)
	{
		for (int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			if(j + 1 < mesh->getNumCells()[1])
			{
				j_plus = j + 1;
			}
			else
			{
				j_plus = j;
			}

			if(j - 1 > -1)
			{
				j_minus = j - 1;
			}
			else
			{
				j_minus = j;
			}

			if(i + 1 < mesh->getNumCells()[0])
			{
				i_plus = i + 1;
			}
			else
			{
				i_plus = i;
			}

			if(i - 1 > -1)
			{
				i_minus = i - 1;
			}
			else
			{
				i_minus = i;
			}

			enstrophy_power += h * h * G * (R[i][j_plus] * sin(k * y[j_plus]) - R[i][j_minus] * sin(k * y[j_minus]))/h * .5 * ((U[i][j_plus] - U[i][j_minus])/h * .5 - (V[i_plus][j] - V[i_minus][j])/h * .5);
		}
	}

	return enstrophy_power;
}
