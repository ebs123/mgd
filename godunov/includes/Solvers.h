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

	//*_with_dummy_reconstr - first index start from 0 (left border) to 4 (bottom border) by clockwise bypass
	void flowDataReconstruction(int* numcells, double** R_with_dummy, double** U_with_dummy, double** V_with_dummy, 
							double** P_with_dummy, double*** R_with_dummy_reconstr, double*** U_with_dummy_reconstr, 
							double*** V_with_dummy_reconstr, double*** P_with_dummy_reconstr, CMeshGenerator* mesh);
	void getFlowDataGradient(double** R_with_dummy, double** U_with_dummy, double** V_with_dummy, 
							double** P_with_dummy, double*** R_with_dummy_grad, double*** U_with_dummy_grad, 
							double*** V_with_dummy_grad, double*** P_with_dummy_grad, CMeshGenerator* mesh);
	void flowVariablesLimiters(double** R_with_dummy, double** U_with_dummy, double** V_with_dummy, double** P_with_dummy, 
							  double*** R_limiter, double*** U_limiter, double*** V_limiter, double*** P_limiter, 
							  CMeshGenerator* mesh);
	void getCellsFlowsSecondOrder(int* numcells, CMeshGenerator* mesh, double** R, double** U, double** V, double** P, double*** dss, double*** uss, 
							 double*** vss, double*** pss);

	~CSolvers(void);
};
