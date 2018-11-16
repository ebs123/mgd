#pragma once
#include "boundary_conds.h"
#include "mesh_generator.h"
#include <algorithm>
#include <sstream>

class CSolvers
{
private:
	double **R;//[i][j] - [x][y]
	double **U;
	double **V;
	double **P;

public:
	CSolvers(int num_cells_x, int num_cells_y);

	void linearSolver(int* numcells, double** R, double** U, double** V, double** P, double*** dss, double*** uss, 
							 double*** vss, double*** pss);
	double getTimeIncrement(double** U, double** V, CMeshGenerator* mesh);
	void solve();

	~CSolvers(void);
};
