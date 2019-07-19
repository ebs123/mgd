#pragma once
#include "boundary_conds.h"

#include "Output.h"
#include "NTracer.h"
#include <algorithm>
#include <sstream>

class CSolvers
{
private:

public:
	CSolvers();

	void linearSolver(int* numcells, double** R, double** U, double** V, double** P, double*** dss, double*** uss, 
							 double*** vss, double*** pss);
	double getTimeIncrement(double** R, double** U, double** V, double** P, CMeshGenerator* mesh);
	void solve(double*** V_init, int n_steps, int n_save, CMeshGenerator* mesh);
	void getGradientsOfFlowData(double*** flow_data_dummy, double*** flow_data);

	~CSolvers(void);
};
