#pragma once
#include "definitions.h"
#include <math.h>

enum boundaryType { SLIP, PERIODIC, OUTFLOW, INFLOW };

class CBoundaryConds
{
private:
	boundaryType x_bound_type, y_bound_type;
	double ***flow_data_dummy;

public:
	CBoundaryConds(void);
	~CBoundaryConds(void);

	void boundary_flows(double **R, double **U, double **V, double **P, double*** dss, double*** uss, double*** vss, 
					double*** pss, int *numcells);
	void linear(double dl, double ul, double vl, double pl, double dr, double ur, double vr, double pr, 
							double &d, double &u, double &v, double &p, char direction);
	double*** getFlowDataWithDummy(double **R, double **U, double **V, double **P, int mesh_size_x, int mesh_size_y); //return array(**R, **U, **V, **P)
};
