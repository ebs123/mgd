#pragma once
#include "NMgdMethods.h"
//#include "NTracer.h"

#ifndef isnan 
#define isnan(x) ((x)!=(x)) 
#endif

class NCalcs
{
public:
	NCalcs(void);
	~NCalcs(void);

	void resid(double ***U, double ***Re,
		const double *hr, const double *hz, double ***Rv);
	void restrict(double ***U, double &t, const double *hr,
		const double *hz);
	void solve(double ***U, double &R2, double &R3, double &R4, double &R5, double &tau, const double *hr,
		const double *hz);

};

