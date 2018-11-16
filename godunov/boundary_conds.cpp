#include "includes\boundary_conds.h"

CBoundaryConds::CBoundaryConds(void)
{
}

CBoundaryConds::~CBoundaryConds(void)
{
}

void CBoundaryConds::boundary_flows(double **R, double **U, double **V, double **P, double*** dss, double*** uss, double*** vss, 
					double*** pss, double *numcells)
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
	for(int j = 0; j < numcells[1]; j++)
	{
		linear(R[numcells[0] - 1][j], U[numcells[0] - 1][j], V[numcells[0] - 1][j], P[numcells[0] - 1][j], R[0][j], U[0][j], 
			V[0][j], P[0][j], dss[0][0][j], uss[0][0][j], vss[0][0][j], pss[0][0][j]);
		linear(R[numcells[0] - 1][j], U[numcells[0] - 1][j], V[numcells[0] - 1][j], P[numcells[0] - 1][j], R[0][j], U[0][j], 
			V[0][j], P[0][j], dss[0][numcells[0]][j], uss[0][numcells[0]][j], vss[0][numcells[0]][j], pss[0][numcells[0]][j]);
	}

	for(int i = 0; i < numcells[0]; i++)
	{
		linear(R[i][numcells[1] - 1], U[i][numcells[1] - 1], V[i][numcells[1] - 1], P[i][numcells[1] - 1], R[i][0], U[i][0], 
			V[i][0], P[i][0], dss[1][i][0], uss[1][i][0], vss[1][i][0], pss[1][i][0]);
		linear(R[i][numcells[1] - 1], U[i][numcells[1] - 1], V[i][numcells[1] - 1], P[i][numcells[1] - 1], R[i][0], U[i][0], 
			V[i][0], P[i][0], dss[1][i][numcells[1]], uss[1][i][numcells[1]], vss[1][i][numcells[1]], pss[1][i][numcells[1]]);
	}
#else
	for(int j = 0; j < numcells[1]; j++)
	{
		linear(R[0][j], U[0][j], V[0][j], P[0][j], R[0][j], U[0][j], V[0][j], P[0][j], dss[0][0][j], uss[0][0][j], vss[0][0][j], pss[0][0][j], 'x');
		linear(R[numcells[0] - 1][j], U[numcells[0] - 1][j], V[numcells[0] - 1][j], P[numcells[0] - 1][j], R[numcells[0] - 1][j], U[numcells[0] - 1][j], 
			V[numcells[0] - 1][j], P[numcells[0] - 1][j], dss[0][numcells[0]][j], uss[0][numcells[0]][j], vss[0][numcells[0]][j], pss[0][numcells[0]][j], 'x');
	}

	for(int i = 0; i < numcells[0]; i++)
	{
		linear(R[i][0], U[i][0], V[i][0], P[i][0], R[i][0], U[i][0], V[i][0], P[i][0], dss[1][i][0], uss[1][i][0], vss[1][i][0], pss[1][i][0], 'y');
		linear(R[i][numcells[1] - 1], U[i][numcells[1] - 1], V[i][numcells[1] - 1], P[i][numcells[1] - 1], R[i][numcells[1] - 1], U[i][numcells[1] - 1], 
			V[i][numcells[1] - 1], P[i][numcells[1] - 1], dss[1][i][numcells[1]], uss[1][i][numcells[1]], vss[1][i][numcells[1]], pss[1][i][numcells[1]], 'y');
	}
#endif
}
//direction == 'x' or 'y'
void CBoundaryConds::linear(double dl, double ul, double vl, double pl, double dr, double ur, double vr, double pr, 
							double &d, double &u, double &v, double &p, char direction) // передача параметров по ссылке
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
		if(direction == 'x')
		{
			bigP = (ul - ur + pl / (dl*cl) + pr / (dr*cr)) / (hl + hr);
			bigU = (dl*cl*ul + dr*cr*ur + pl - pr) / (dl*cl + dr*cr);
			R3 = dl - dl / cl * (bigU - ul);
			R4 = dr + dr / cr * (bigU - ur);
			if (bigU > 0)
			{
				bigR = R3;
				bigV = vl;
			}
			else
			{
				bigR = R4;
				bigV = vr;
			}
		}
		else
		{
			bigP = (vl - vr + pl / (dl*cl) + pr / (dr*cr)) / (hl + hr);
			bigV = (dl*cl*vl + dr*cr*vr + pl - pr) / (dl*cl + dr*cr);
			R3 = dl - dl / cl * (bigV - vl);
			R4 = dr + dr / cr * (bigV - vr);
			if (bigV > 0)
			{
				bigR = R3;
				bigU = ul;
			}
			else
			{
				bigR = R4;
				bigU = ur;
			}
		}
		/*	if (bigU >= 0) bigS = pl / pow(dl, GAMMA);
		else bigS = pr / pow(dr, GAMMA);
		help = bigP / bigS;
		bigR = pow(help, 1.0 / GAMMA);*/
	}

	u = bigU;
	p = bigP;
	d = bigR;
	v = bigV;

}


