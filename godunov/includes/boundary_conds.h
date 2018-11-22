#pragma once
#include "definitions.h"
#include <math.h>

class CBoundaryConds
{
public:
	CBoundaryConds(void);
	~CBoundaryConds(void);

	void boundary_flows(double **R, double **U, double **V, double **P, double*** dss, double*** uss, double*** vss, 
					double*** pss, int *numcells);

	void linear(double dl, double ul, double vl, double pl, double dr, double ur, double vr, double pr, 
							double &d, double &u, double &v, double &p, char direction);
};
