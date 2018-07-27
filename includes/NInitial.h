#pragma once
#ifndef NINITIAL_H
#define NINITIAL_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <math.h>
#include "methods.h"
#include <map>
#include <iostream>
#include <vector>

using namespace std;

class NInitial
{
private:
	static int Nsave, Nstep, Index, RZONE, Ncomp, xmax, ymax;
	static double Cp, kappa, RFAKTOR, roin, uphiin, uzin, urin, pin, sigma, deltaz, deltar, R1, mu, H[3];
	static multimap<string, string> xmlTagOptionsMap;
	static map<string, double> initialData;
	static vector<string> tagsList;
	static const double pi;
	static const char* slip_name;
	static const char* periodic_name;

protected:
	static int parseXMLTag(string &str);

public:
	NInitial(void);
	~NInitial(void);

	double domega(const double &r);
	double ad(const double &r);
	double adr(const double &r, const double &z);
	double U_phi(const double &z, const double &r, double &omega);
	void initial(double ***U, const double *hr, const double *hz);
	static void loadInitialData(const string &patch);
	static int getNsave();
	static int getNstep();
	static int getIndex();
	static int getRZONE();
	static int getNcomp();
	static int get_xmax();
	static int get_ymax();
	static double getCp();
	static double get_kappa();
	static double getRFAKTOR();
	static double get_roin();
	static double get_uphiin();
	static double get_uzin();
	static double get_urin();
	static double get_pin();
	static double get_sigma();
	static double get_deltaz();
	static double get_deltar();
	static double getR1();
	static double get_mu();
	static double get_pi();
	static double* getH();

	static char* getPeriodicName();
	static char* getSlipName();
};

#endif NINITIAL_H