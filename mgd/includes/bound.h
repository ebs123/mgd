#ifndef BOUND_H
#define BOUND_H


#include <iostream>	
#include <math.h>
#include "NInitial.h"
//#include "methods.h"
//#define Ncomp 5
//#define pin 13
//#define xmax 300
//#define ymax 600
//#define sigma 0.4
//#define deltar 0.5
//#define deltaz 1.0
//#define R1    .5
//#define mu 1.811

using namespace std;

class NBoundaryCond
{
private:
public:
	virtual ~NBoundaryCond();
	virtual void slip(const double *U, double *U1, char *dir);
	virtual void slipinner(const double *U, double *U1);
	virtual void periodic(const int &i, const int &j, double ***U, double *U1, char *dir);
	void inflow(const int &i, const double &p, const double &ro, const double &ur, double uz, const double &uphi,
			double *U1, const double *hr);
	void outflow(const double *U, double *U1);
	void slipVisc(const double *U, double *U1, char *dir);
	void slipinnerVisc(const double *U, double *U1);
};

class NMgdBoundaryCond :
	public NBoundaryCond
{
public:

	NMgdBoundaryCond(void);
	~NMgdBoundaryCond(void);

	virtual void slip(const double *U, double *U1, char *dir);
	virtual void slipinner(const double *U, double *U1);
	virtual void periodic(const int &i, const int &j, double ***U, double *U1, char *dir);
};


#endif //BOUND_H