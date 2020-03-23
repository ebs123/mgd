#pragma once
#include "methods.h"
#include "NTracer.h"


class NMgdMethods :
	public NMethods
{
protected:
	vector<double> getSVariables(const double *Ur, const double *Ul);
	vector<vector<double> > roeMatrix(const double *Ur, const double *Ul, const char *dir, const char *dirn);

public:
	virtual void Frp(const double *U, double *F);
	virtual void Fz(const double *U, double *F);
	double* f_cusp(const double *U, const string &direct);
	double* P_cusp(const double *U, const string &direct);
	double* psi_cusp(const double *U);
	virtual void Roe(const double *Ur, const double *Ul, double *Flux, const char *dir, const char *dirn);
	virtual void Sourse(const double *U, double *S);
	double energy(const double *U);
	double*** korrectH(const double ***U, const double *hr, const double *hz, double &tau);

	NMgdMethods(void);
	~NMgdMethods(void);
};
