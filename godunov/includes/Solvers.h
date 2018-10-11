#pragma once

class CSolvers
{
private:
	double **R;
	double **U;
	double **V;
	double **P;

public:
	CSolvers(int num_cells_x, int num_cells_y);

	void linearSolver(int numcells, double* R, double* U, double* P, double* dss, double* uss, double* pss);
	void solve();

	~CSolvers(void);
};
