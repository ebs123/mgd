#include "includes\boundary_conds.h"

CBoundaryConds::CBoundaryConds(boundaryType _x_bound_type, boundaryType _y_bound_type) : x_bound_type(_x_bound_type), 
	y_bound_type(_y_bound_type)
{
	flow_data_dummy = 0;
}

CBoundaryConds::~CBoundaryConds(void)
{
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < m_mesh_size_x + 4; j++)
		{
			delete[] flow_data_dummy[i][j];
		}

	for(int i = 0; i < 4; i++)
	{
		delete[] flow_data_dummy[i];
	}

	delete[] flow_data_dummy;
}

void CBoundaryConds::boundary_flows(double **R, double **U, double **V, double **P, double*** dss, double*** uss, double*** vss, 
					double*** pss, int *numcells)
{
	double c0, s0, cn, sn;

	for(int j = 0; j < numcells[1]; j++)
	{
		linear(R[0][j], -U[0][j], V[0][j], P[0][j], R[0][j], U[0][j], V[0][j], P[0][j], dss[0][0][j], uss[0][0][j], vss[0][0][j], pss[0][0][j], 'x');
		linear(R[numcells[0] - 1][j], U[numcells[0] - 1][j], V[numcells[0] - 1][j], P[numcells[0] - 1][j], R[numcells[0] - 1][j], -U[numcells[0] - 1][j], 
			V[numcells[0] - 1][j], P[numcells[0] - 1][j], dss[0][numcells[0]][j], uss[0][numcells[0]][j], vss[0][numcells[0]][j], pss[0][numcells[0]][j], 'x');
	}

	for(int i = 0; i < numcells[0]; i++)
	{
		linear(R[i][0], U[i][0], -V[i][0], P[i][0], R[i][0], U[i][0], V[i][0], P[i][0], dss[1][i][0], uss[1][i][0], vss[1][i][0], pss[1][i][0], 'y');
		linear(R[i][numcells[1] - 1], U[i][numcells[1] - 1], V[i][numcells[1] - 1], P[i][numcells[1] - 1], R[i][numcells[1] - 1], U[i][numcells[1] - 1], 
			-V[i][numcells[1] - 1], P[i][numcells[1] - 1], dss[1][i][numcells[1]], uss[1][i][numcells[1]], vss[1][i][numcells[1]], pss[1][i][numcells[1]], 'y');
	}
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

double*** CBoundaryConds::getFlowDataWithDummy(double **R, double **U, double **V, double **P, int mesh_size_x, int mesh_size_y)
{
	if(flow_data_dummy != 0)
	{
		return flow_data_dummy;
	}

	m_mesh_size_x = mesh_size_x;

	flow_data_dummy = new double**[4];
	for(int i = 0; i < 4; i++)
	{
		flow_data_dummy[i] = new double*[mesh_size_x + 4];
	}
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < mesh_size_x + 4; j++)
		{
			flow_data_dummy[i][j] = new double[mesh_size_y + 4];
		}

	for(int i = 0; i < mesh_size_x; i++)
		for(int j = 0; j < mesh_size_y; j++)
		{
			flow_data_dummy[0][i + 2][j + 2] = R[i][j];
			flow_data_dummy[1][i + 2][j + 2] = U[i][j];
			flow_data_dummy[2][i + 2][j + 2] = V[i][j];
			flow_data_dummy[3][i + 2][j + 2] = P[i][j];
		}

	for(int j = 0; j < mesh_size_y; j++)
	{
		if(x_bound_type == OUTFLOW)
		{
			flow_data_dummy[0][0][j + 2] = R[0][j];
			flow_data_dummy[1][0][j + 2] = U[0][j];
			flow_data_dummy[2][0][j + 2] = V[0][j];
			flow_data_dummy[3][0][j + 2] = P[0][j];

			flow_data_dummy[0][1][j + 2] = R[0][j];
			flow_data_dummy[1][1][j + 2] = U[0][j];
			flow_data_dummy[2][1][j + 2] = V[0][j];
			flow_data_dummy[3][1][j + 2] = P[0][j];

			flow_data_dummy[0][mesh_size_x + 2][j + 2] = R[mesh_size_x - 1][j];
			flow_data_dummy[1][mesh_size_x + 2][j + 2] = U[mesh_size_x - 1][j];
			flow_data_dummy[2][mesh_size_x + 2][j + 2] = V[mesh_size_x - 1][j];
			flow_data_dummy[3][mesh_size_x + 2][j + 2] = P[mesh_size_x - 1][j];

			flow_data_dummy[0][mesh_size_x + 3][j + 2] = R[mesh_size_x - 1][j];
			flow_data_dummy[1][mesh_size_x + 3][j + 2] = U[mesh_size_x - 1][j];
			flow_data_dummy[2][mesh_size_x + 3][j + 2] = V[mesh_size_x - 1][j];
			flow_data_dummy[3][mesh_size_x + 3][j + 2] = P[mesh_size_x - 1][j];
		}
		else if(x_bound_type == PERIODIC)
		{
			flow_data_dummy[0][0][j + 2] = R[1][j];
			flow_data_dummy[1][0][j + 2] = U[1][j];
			flow_data_dummy[2][0][j + 2] = V[1][j];
			flow_data_dummy[3][0][j + 2] = P[1][j];

			flow_data_dummy[0][1][j + 2] = R[0][j];
			flow_data_dummy[1][1][j + 2] = U[0][j];
			flow_data_dummy[2][1][j + 2] = V[0][j];
			flow_data_dummy[3][1][j + 2] = P[0][j];

			flow_data_dummy[0][mesh_size_x + 2][j + 2] = R[mesh_size_x - 1][j];
			flow_data_dummy[1][mesh_size_x + 2][j + 2] = U[mesh_size_x - 1][j];
			flow_data_dummy[2][mesh_size_x + 2][j + 2] = V[mesh_size_x - 1][j];
			flow_data_dummy[3][mesh_size_x + 2][j + 2] = P[mesh_size_x - 1][j];

			flow_data_dummy[0][mesh_size_x + 3][j + 2] = R[mesh_size_x - 2][j];
			flow_data_dummy[1][mesh_size_x + 3][j + 2] = U[mesh_size_x - 2][j];
			flow_data_dummy[2][mesh_size_x + 3][j + 2] = V[mesh_size_x - 2][j];
			flow_data_dummy[3][mesh_size_x + 3][j + 2] = P[mesh_size_x - 2][j];
		}
	}

	for(int i = 0; i < mesh_size_x; i++)
	{
		if(y_bound_type == OUTFLOW)
		{
			flow_data_dummy[0][i + 2][0] = R[i][0];
			flow_data_dummy[1][i + 2][0] = U[i][0];
			flow_data_dummy[2][i + 2][0] = V[i][0];
			flow_data_dummy[3][i + 2][0] = P[i][0];

			flow_data_dummy[0][i + 2][1] = R[i][0];
			flow_data_dummy[1][i + 2][1] = U[i][0];
			flow_data_dummy[2][i + 2][1] = V[i][0];
			flow_data_dummy[3][i + 2][1] = P[i][0];

			flow_data_dummy[0][i + 2][mesh_size_y + 2] = R[i][mesh_size_y - 1];
			flow_data_dummy[1][i + 2][mesh_size_y + 2] = U[i][mesh_size_y - 1];
			flow_data_dummy[2][i + 2][mesh_size_y + 2] = V[i][mesh_size_y - 1];
			flow_data_dummy[3][i + 2][mesh_size_y + 2] = P[i][mesh_size_y - 1];

			flow_data_dummy[0][i + 2][mesh_size_y + 3] = R[i][mesh_size_y - 1];
			flow_data_dummy[1][i + 2][mesh_size_y + 3] = U[i][mesh_size_y - 1];
			flow_data_dummy[2][i + 2][mesh_size_y + 3] = V[i][mesh_size_y - 1];
			flow_data_dummy[3][i + 2][mesh_size_y + 3] = P[i][mesh_size_y - 1];
		}
		else if(y_bound_type == PERIODIC)
		{
			flow_data_dummy[0][i + 2][0] = R[i][1];
			flow_data_dummy[1][i + 2][0] = U[i][1];
			flow_data_dummy[2][i + 2][0] = V[i][1];
			flow_data_dummy[3][i + 2][0] = P[i][1];

			flow_data_dummy[0][i + 2][1] = R[i][0];
			flow_data_dummy[1][i + 2][1] = U[i][0];
			flow_data_dummy[2][i + 2][1] = V[i][0];
			flow_data_dummy[3][i + 2][1] = P[i][0];

			flow_data_dummy[0][i + 2][mesh_size_y + 2] = R[i][mesh_size_y - 1];
			flow_data_dummy[1][i + 2][mesh_size_y + 2] = U[i][mesh_size_y - 1];
			flow_data_dummy[2][i + 2][mesh_size_y + 2] = V[i][mesh_size_y - 1];
			flow_data_dummy[3][i + 2][mesh_size_y + 2] = P[i][mesh_size_y - 1];

			flow_data_dummy[0][i + 2][mesh_size_y + 3] = R[i][mesh_size_y - 2];
			flow_data_dummy[1][i + 2][mesh_size_y + 3] = U[i][mesh_size_y - 2];
			flow_data_dummy[2][i + 2][mesh_size_y + 3] = V[i][mesh_size_y - 2];
			flow_data_dummy[3][i + 2][mesh_size_y + 3] = P[i][mesh_size_y - 2];
		}
	}

	return flow_data_dummy;
}
