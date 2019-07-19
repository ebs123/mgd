#pragma once
#include "mesh_generator.h"
#include "definitions.h"
#include <math.h>

class CFlowParameters
{
public:
	CFlowParameters(void);
	double getKineticEnergySpectrum(double k_x, double k_y, double** U, double** V, CMeshGenerator* mesh);

public:
	~CFlowParameters(void);
};
