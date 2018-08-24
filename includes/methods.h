#pragma once
#ifndef METHODS_H
#define METHODS_H

#include "NMgdBoundaryCond.h"
#include <vector>

#include "NArrPacker.h"
#define DUMMY_NUM 1
#define DUMMY_NUM2 3

class NMethods
{
protected:
	double sign(const double &val);

public:
	virtual ~NMethods();
	void LoadGrid(double *hr, double *hz);
	double zij1(const int &j, const double &h0);
	double rij1(const int &i, const double *hr);//h0 - шаг первой €чейки слева
	double press(const double *U);
	double sound(const double *U);
	double sound(const double *U, const char *dir);
	double MakhNumber(const double *U);
	double entalpy(const double *U);
	void get_charact_variables(const double *U, double *V);
	virtual void Frp(const double *U, double *F);
	virtual void Fz(const double *U, double *F);
	virtual void Roe(const double *Ur, const double *Ul, double *Flux, const char *dir, const char *dirn);
	virtual void Sourse(const double *U, double *S);
	void get_grad1(const int i, const int j, double ***U, double **grad, const char *bz, const char *br, const double *hr, const double *hz);
	void get_grad1Visc(const int i, const int j, const double ***U, double **grad, const char *bz, const char *br,
		const double *hr, const double *hz);
	void get_grad2(const int i, const int j, const double ***U, double **grad, const char *bz, const char *br,
		const double *hr, const double *hz);
	void getGradH(const int i, const int j, const double ***U, const double *hr,
		const double *hz, const char *bz, const char *br, double gradH[2]);
	void SourseVisc(const double *U, double *S, double (&tau)[3][3], const double Hr);
	void stress(const int i, const int j, const double ***U, double (&tau)[3][3], char *bz, char *br,
		const double *hr, const double *hz);
	void Fzvisc(const double *U, double *F, double tau[3][3], const double Hz);
	void Frvisc(const double *U, double *F, double tau[3][3], const double Hr);
	void limi(const int i, const int j, const double ***U, double *psi, char *dir, char *bz, char *br,
		const double *hr, const double *hz, const char *dirn);
	void limibound(const int i, const int j, const double ***U, double *psi, char *dir, char *bz, char *br,
		const double *hr, const double *hz, const char *dirn);
	void limiter(const int i, const int j, double ***U, double *psi, char *bz, char *br,
		const double *hr, const double *hz, const char *dirn);

};


#endif //METHODS_H