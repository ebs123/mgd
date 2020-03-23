#pragma once
#ifndef BOUND_CPP
#define BOUND_CPP

#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\bound.h"
#else
#include "includes/bound.h"
#endif

NBoundaryCond::~NBoundaryCond()
{
	;
}
void NBoundaryCond::slip(const double *U, double *U1, char *dir)
{
	if (dir == "r")
	{
		U1[1] = -U[1];
		U1[2] = U[2];
	}
	else if(dir == "z")
	{
		U1[1] = U[1];
		U1[2] = -U[2];
	}
	else
	{
		cout << "err in slip dir";
		exit(1);
	}

	U1[0] = U[0];
	U1[3] = U[3];
	U1[4] = U[4];
}
void NBoundaryCond::slipinner(const double *U, double *U1)
{
	U1[0] = U[0];
	U1[1] = -U[1];
	U1[2] = U[2];
	U1[3] = U[3];
	U1[4] = U[4];
}
void NBoundaryCond::periodic(const int &i, const int &j, double ***U, double *U1, char *dir)
{
	if(dir == "z")
	{ 
		const double p1 = 10;

		if(j < 10)
		{
			//U1[Ncomp-1]=p1;
			for(int k = 0; k < NInitial::getNcomp(); k++)
				U1[k] = U[i][NInitial::get_ymax() - j - 1][k];
		}
		else if(j > NInitial::get_ymax() - 11)
		{
			//U1[Ncomp-1]=pin;
			for(int k = 0; k < NInitial::getNcomp(); k++)
				U1[k] = U[i][NInitial::get_ymax() - 1 - j][k];
		}
		else
		{
			char *a;
			cout << j << ' ' << "err in bound z" << '\n';
			cin >> a;
			exit(1);
		}

	}

	if(dir == "r")
	{
		if(i < 100)
			for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				U1[k] = U[NInitial::get_xmax() - i - 1][j][k];
		else if(i > NInitial::get_xmax() - 11)
			for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				U1[k] = U[NInitial::get_xmax() - 1 - i][j][k];
		else
		{
			char *a;
			cout << i << ' ' << "err in bound r" << '\n';
			cin >> a;
			exit(1);
		}

	}
	
}
void NBoundaryCond::inflow(const int &i, const double &p, const double &ro, const double &ur, double uz, const double &uphi,
			double *U1, const double *hr)
{
	double p0 = 10, p1 = 13, rij, R2;
	NMethods methods;

	rij = methods.rij1(i, hr);
	R2 = NInitial::getR1() + NInitial::get_deltar();

	uz = 0.25 * (p1 - p0)/(NInitial::get_mu() * ro * NInitial::get_deltaz()) * (pow(rij, 2) - pow(NInitial::getR1(), 2) - (pow(R2, 2) - pow(NInitial::getR1(), 2)) * 
		log(rij/NInitial::getR1())/log(R2/NInitial::getR1()))/*ad(rij)*/;

	U1[0] = ro;
	U1[1] = ro * ur;
	U1[2] = ro * uz;
	U1[3] = ro * uphi;
	U1[4] = ro * (p/(NInitial::get_sigma() * ro) + .5 * (pow(ur, 2) + pow(uz, 2) + pow(uphi, 2)));
}
void NBoundaryCond::outflow(const double *U, double *U1)
{

	for(int i = 0; i < NInitial::getNcomp(); i++)
	U1[i] = U[i];
}
void NBoundaryCond::slipVisc(const double *U, double *U1, char *dir)
{
	U1[0] = U[0];
	//U1[1]=-U[1];
	//U1[2]=-U[2];

	if (dir == "r")
	{
		U1[1] = -U[1];
		U1[2] = -U[2];
		U1[3] = .5 * (NInitial::getR1() + NInitial::get_deltar());//1 == omega, in initial conditions
	}
	else if(dir == "z")
	{
		U1[1] = -U[1];
		U1[2] = -U[2];
		U1[3] = -U[3];
	}
	else
	{
		cout << "err in slip dir";
		exit(1);
	}

	//U1[3]=U[3];
	//U1[3]=U[3];
	U1[4] = U[4] - .5 * (pow(U[1], 2) + pow(U[2], 2) + pow(U[3], 2))/U[0] + 
		.5 * (pow(U1[1], 2) + pow(U1[2], 2) + pow(U1[3], 2))/U1[0];
}
void NBoundaryCond::slipinnerVisc(const double *U, double *U1)
{
	U1[0] = U[0];
	U1[1] = -U[1];
	U1[2] = -U[2];
	U1[3] = .5 * NInitial::getR1();//1 == omega, in initial conditions
	U1[4] = U[4] - .5*(pow(U[1], 2) + pow(U[2], 2) + pow(U[3], 2))/U[0] + 
		.5*(pow(U1[1], 2) + pow(U1[2], 2) + pow(U1[3], 2))/U1[0];
}

#endif //BOUND_CPP