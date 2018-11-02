#include "boundary_conds.h"

CBoundaryConds::CBoundaryConds(void)
{
}

CBoundaryConds::~CBoundaryConds(void)
{
}

void CBoundaryConds::boundary_flows(double **R, double **U, double **V, double **P, double*** dss, double*** uss, double*** vss, 
					double*** pss)
{
	double c0, s0, cn, sn;

#if (PROBLEM==18 || PROBLEM==20)
	// set pressure
	pss[0] = initial_pressure(0.0);
	//	pss[numcells] = initial_pressure(LENGTH);

	c0 = sqrt(GAMMA*P[0] / R[0]);
	//	cn = sqrt(GAMMA*P[numcells] / R[numcells]);

	double l0_const = U[0] - P[0] / (R[0] * c0);
	//double rn_const = U[numcells] + P[numcells] / (R[numcells] * cn);

	uss[0] = l0_const + pss[0] / (R[0] * c0);
	//	uss[numcells] = rn_const - pss[numcells] / (R[numcells] * cn);

	s0 = log(P[0] / pow(R[0], GAMMA));
	//	sn = log(P[numcells] / pow(R[numcells], GAMMA));

	//dss[0] = pow(pss[0] / s0, 1.0 / GAMMA);
	//	dss[numcells] = pow(pss[numcells] / sn, 1.0 / GAMMA);

	//R3 = dl - dl / cl * (bigU - ul);
	dss[0] = R[0] + R[0] / c0 * (uss[0] - U[0]);

	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);

#elif (PROBLEM == 19)

	uss[0] = timer;
	c0 = sqrt(GAMMA*P[0] / R[0]);

	double l0_const = U[0] - P[0] / (R[0] * c0);
	double rn_const = U[numcells] + P[numcells] / (R[numcells] * Cn);

	pss[0] = (uss[0] - l0_const)*(R[0] * c0);

	s0 = log(P[0] / pow(R[0], GAMMA));
	dss[0] = pow(pss[0] / s0, 1.0 / GAMMA);


#elif (PROBLEM == 3)
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[0], U[0], P[0], dss[numcells], uss[numcells], pss[numcells]);
#else
	linear(R[0], U[0], P[0], R[0], U[0], P[0], dss[0], uss[0], pss[0]);
	linear(R[numcells - 1], U[numcells - 1], P[numcells - 1], R[numcells - 1], U[numcells - 1], P[numcells - 1], dss[numcells], uss[numcells], pss[numcells]);
#endif
}

void CBoundaryConds::linear(double dl, double ul, double vl, double pl, double dr, double ur, double vr, double pr, 
							double &d, double &u, double &v, double &p) // передача параметров по ссылке
{
	double bigU, bigV, bigP, bigS, bigR, cl, cr, help, hl, hr, R3, R4;

	cl = sqrt(GAMMA*pl / dl);
	cr = sqrt(GAMMA*pr / dr);
	hl = 1.0 / (dl*cl);
	hr = 1.0 / (dr*cr);

	if (ul > cl)
	{
		bigP = pl;
		bigU = ul;
		bigR = dl;
		bigV = vl;
	}
	else if (ur < -cr)
	{
		bigP = pr;
		bigU = ur;
		bigR = dr;
		bigV = vr;
	}
	else
	{
		bigP = (ul - ur + pl / (dl*cl) + pr / (dr*cr)) / (hl + hr);
		bigU = (dl*cl*ul + dr*cr*ur + pl - pr) / (dl*cl + dr*cr);
		/*	if (bigU >= 0) bigS = pl / pow(dl, GAMMA);
		else bigS = pr / pow(dr, GAMMA);
		help = bigP / bigS;
		bigR = pow(help, 1.0 / GAMMA);*/
		R3 = dl - dl / cl * (bigU - ul);
		R4 = dr + dr / cr * (bigU - ur);
		if (bigU > 0) bigR = R3;
		else bigR = R4;
	}

	u = bigU;
	p = bigP;
	d = bigR;

}


