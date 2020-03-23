//#include "StdAfx.h"
#include "includes\nmgdboundarycond.h"

NMgdBoundaryCond::NMgdBoundaryCond(void)
{
}

NMgdBoundaryCond::~NMgdBoundaryCond(void)
{
}

void NMgdBoundaryCond::slip(const double *U, double *U1, char *dir)
{
	NBoundaryCond::slip(U, U1, dir);

	if (dir == "r")
	{
		double *H = NInitial::getH();
		U1[5] = U[5];
		U1[6] = U[6];
		U1[7] = U[7];
	}
	else if(dir == "z")
	{
//		double *H = NInitial::getH();
		U1[5] = U[5];
		U1[6] = U[6];
		U1[7] = U[7];
	}
	else
	{
		cout << "err in slip dir";
		exit(1);
	}
}
void NMgdBoundaryCond::slipinner(const double *U, double *U1)
{
	NBoundaryCond::slipinner(U, U1);

	//double *H = NInitial::getH();
	U1[5] = U[5];
	U1[6] = U[6];
	U1[7] = U[7];
}
void NMgdBoundaryCond::periodic(const int &i, const int &j, double ***U, double *U1, char *dir)
{
	NBoundaryCond::periodic(i, j, U, U1, dir);
}
