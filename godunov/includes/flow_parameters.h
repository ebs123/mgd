#pragma once
#include "mesh_generator.h"
#include "definitions.h"
#include <math.h>

class CFlowParameters
{
private:
	double *E;
	double **S;
	double **S_last_step;
	double **S_diff;

public:
	int k_x_max;
	int k_y_max;
	CMeshGenerator *mesh;

public:
	CFlowParameters(CMeshGenerator *grid);
	double* getKineticEnergySpectrum(double** U, double** V);
	double** getFlowEntropy(double** R, double** P);
	double** getFlowEntropyDiff();
	double getEnergyPower(double** U, double** R, const double G, const double k);
	double getEnstrophyPower(double** U, double** V, double** R, const double G, const double k);

	int k_max();
	void saveFlowEntropyCurrTimeStep(double** R, double** P);

public:
	~CFlowParameters(void);
};

inline int CFlowParameters::k_max()
{
	return (int)sqrt((double)(k_x_max * k_x_max + k_y_max * k_y_max));
}

