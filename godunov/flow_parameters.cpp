#include "includes\flow_parameters.h"

CFlowParameters::CFlowParameters(void)
{
}

CFlowParameters::~CFlowParameters(void)
{
}

double CFlowParameters::getKineticEnergySpectrum(double k_x, double k_y, double** U, double** V, CMeshGenerator* mesh)
{
	double e_1, e_2, e_3, e_4;
	double *x, *y;
	double dx, dy;

	x = mesh->getMeshComponent(0);
	y = mesh->getMeshComponent(1);
	dx = mesh->getMeshStep(0);
	dy = mesh->getMeshStep(1);

	e_1 = 0;
	e_2 = 0;
	e_3 = 0;
	e_4 = 0;
	for(int i = 0; i < mesh->getNumCells()[0]; i++)
	{
		if(x[i] > 2 * pi)
		{
			break;
		}

		for(int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			if(y[j] > 2 * pi)
			{
				break;
			}

			e_1 += 1/pi * (U[i][j] * U[i][j] + V[i][j] * V[i][j]) * .5 * cos(k_x * x[i]) * cos(k_y * y[j]) * dx * dy;
			e_2 += 1/pi * (U[i][j] * U[i][j] + V[i][j] * V[i][j]) * .5 * cos(k_x * x[i]) * sin(k_y * y[j]) * dx * dy;
			e_3 += 1/pi * (U[i][j] * U[i][j] + V[i][j] * V[i][j]) * .5 * sin(k_x * x[i]) * cos(k_y * y[j]) * dx * dy;
			e_4 += 1/pi * (U[i][j] * U[i][j] + V[i][j] * V[i][j]) * .5 * sin(k_x * x[i]) * sin(k_y * y[j]) * dx * dy;
		}
	}

	return sqrt(e_1 * e_1 + e_2 * e_2 + e_3 * e_3 + e_4 * e_4);
}