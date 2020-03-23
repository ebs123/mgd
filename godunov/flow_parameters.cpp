#include "includes\flow_parameters.h"

CFlowParameters::CFlowParameters(void)
{
	k_x_max = 200;
	k_y_max = 200;
	E = new double[(int)sqrt((double)(k_x_max * k_x_max + k_y_max * k_y_max))];
}

CFlowParameters::~CFlowParameters(void)
{
	delete E;
}

double* CFlowParameters::getKineticEnergySpectrum(double** U, double** V)
{
	double U_fourier_i[4], V_fourier_i[4];
	const int N = 398;
	double h, E_kx_ky;
	const int k_max = sqrt((double)(k_x_max * k_x_max + k_y_max * k_y_max));

	for(int i = 0; i < (int)sqrt((double)(k_x_max * k_x_max + k_y_max * k_y_max)); i++)
	{
		E[i] = 0;
	}

	for(int i = 0; i < 4; i++)
	{
		U_fourier_i[i] = 0;
		V_fourier_i[i] = 0;
	}

	h = 2 * pi/(N + 1);

	for(int k_loop = 1; k_loop < k_max + 1; k_loop++)
	{
		for(int k_x = 1; k_x < k_x_max + 1; k_x++)
		{
			for(int k_y = 1; k_y < k_y_max + 1; k_y++)
			{
				double k_mod = sqrt((double)(k_x_max * k_x_max + k_y_max * k_y_max));
				for(int i = 1; i < N + 2; i++)
				{
					for(int j = 1; j < N + 2; j++)
					{
						U_fourier_i[0] += 1/pi * U[i][j] * cos(k_x * h * i * .5) * cos(k_y * h * j * .5) * h * h;
						U_fourier_i[1] += 1/pi * U[i][j] * cos(k_x * h * i * .5) * sin(k_y * h * j * .5) * h * h;
						U_fourier_i[2] += 1/pi * U[i][j] * sin(k_x * h * i * .5) * cos(k_y * h * j * .5) * h * h;
						U_fourier_i[3] += 1/pi * U[i][j] * sin(k_x * h * i * .5) * sin(k_y * h * j * .5) * h * h;

						V_fourier_i[0] += 1/pi * U[i][j] * cos(k_x * h * i * .5) * cos(k_y * h * j * .5) * h * h;
						V_fourier_i[1] += 1/pi * U[i][j] * cos(k_x * h * i * .5) * sin(k_y * h * j * .5) * h * h;
						V_fourier_i[2] += 1/pi * U[i][j] * sin(k_x * h * i * .5) * cos(k_y * h * j * .5) * h * h;
						V_fourier_i[3] += 1/pi * U[i][j] * sin(k_x * h * i * .5) * sin(k_y * h * j * .5) * h * h;
					}
				}

				E_kx_ky = (U_fourier_i[0] * U_fourier_i[0] + U_fourier_i[1] * U_fourier_i[1] + U_fourier_i[2] * U_fourier_i[2] + 
					U_fourier_i[3] * U_fourier_i[3] + V_fourier_i[0] * V_fourier_i[0] + V_fourier_i[1] * V_fourier_i[1] + 
					V_fourier_i[2] * V_fourier_i[2] + V_fourier_i[3] * V_fourier_i[3]) * .5;

				if(abs(k_loop - k_mod) < .5)
				{
					E[k_loop] += E_kx_ky;
				}
			}
		}
	}

	return E;
}