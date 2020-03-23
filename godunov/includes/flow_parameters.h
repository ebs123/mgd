#pragma once
#include "mesh_generator.h"
#include "definitions.h"
#include <math.h>

class CFlowParameters
{
private:
	double *E;
public:
	int k_x_max;
	int k_y_max;

public:
	CFlowParameters(void);
	double* getKineticEnergySpectrum(double** U, double** V);

public:
	~CFlowParameters(void);
};
