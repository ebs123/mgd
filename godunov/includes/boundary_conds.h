#pragma once

class CBoundaryConds
{
public:
	CBoundaryConds(void);
public:
	~CBoundaryConds(void);

private:

	void boundary_flows(double **R, double **U, double **V, double **P, double*** dss, double*** uss, double*** vss, 
					double*** pss, double *numcells);

	void linear(double dl, double ul, double pl, double dr, double ur, double pr, double &d, double &u, double &p);
};
