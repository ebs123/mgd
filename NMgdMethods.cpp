//#include "StdAfx.h"
#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\nmgdmethods.h"
#else
#include "includes/nmgdmethods.h"
#endif


NMgdMethods::NMgdMethods(void)
{
}

NMgdMethods::~NMgdMethods(void)
{
}

void NMgdMethods::Frp(const double *U, double *F)
{
	NMethods::Frp(U, F);
	F[1] += .5 * (pow(U[6], 2) + pow(U[7], 2) - pow(U[5], 2));
	F[2] -= U[5] * U[6];
	F[3] -= U[5] * U[7];
	F[4] += U[1]/U[0] * .5 * (pow(U[5], 2) + pow(U[6], 2) + pow(U[7], 2)) - 
		U[5] * (U[5] * U[1]/U[0] + U[6] * U[2]/U[0] + U[7] * U[3]/U[0]);
	F[5] = 0;
	F[6] = U[1]/U[0] * U[6] - U[2]/U[0] * U[5];
	F[7] = U[1]/U[0] * U[7] - U[3]/U[0] * U[5];
}

void NMgdMethods::Fz(const double *U, double *F)
{
	NMethods::Fz(U, F);
	F[1] -= U[5] * U[6];
	F[2] += .5 * (pow(U[5], 2) + pow(U[7], 2) - pow(U[6], 2));
	F[3] -= U[6] * U[7];
	F[4] += U[2]/U[0] * .5 * (pow(U[5], 2) + pow(U[6], 2) + pow(U[7], 2)) - 
		U[6] * (U[5] * U[1]/U[0] + U[6] * U[2]/U[0] + U[7] * U[3]/U[0]);
	F[5] = U[2]/U[0] * U[5] - U[1]/U[0] * U[6];
	F[6] = 0.;
	F[7] = U[2]/U[0] * U[7] - U[3]/U[0] * U[6];
}

double* NMgdMethods::f_cusp(const double *U, const string &comp)//выделяет память под return массив!
{
	double *ret = new double[NInitial::getNcomp()];

	for(size_t i = 0; i < NInitial::getNcomp(); i++)
	{
		if(i == 5 & !strcmp(comp.c_str(), "r"))
		{
			ret[i] = 0.;
			continue;
		}
		else if(i == 6 & !strcmp(comp.c_str(), "z"))
		{
			ret[i] = 0.;
			continue;
		}
		ret[i] = U[i];
	}

	return ret;
}

double* NMgdMethods::P_cusp(const double *U, const string &comp)//выделяет память под return массив!
{
	double *ret = new double[NInitial::getNcomp()];
	NMethods methods;

	if(!strcmp(comp.c_str(), "r"))
	{
		ret[0] = 0.;
		ret[1] = methods.press(U) + .5 * (pow(U[6], 2) - pow(U[5], 2) + pow(U[7], 2));
		ret[2] = - U[5] * U[6];
		ret[3] = - U[5] * U[7];
		ret[4] = - U[5] * (U[5] * U[1]/U[0] + U[6] * U[2]/U[0] + U[7] * U[3]/U[0]);
		ret[5] = 0.;
		ret[6] = - U[2]/U[0] * U[5];
		ret[7] = - U[3]/U[0] * U[5];
	}
	else if(!strcmp(comp.c_str(), "z"))
	{
		ret[0] = 0.;
		ret[1] = - U[5] * U[6];
		ret[2] = methods.press(U) + .5 * (pow(U[5], 2) - pow(U[6], 2) + pow(U[7], 2));
		ret[3] = - U[6] * U[7];
		ret[4] = - U[6] * (U[5] * U[1]/U[0] + U[6] * U[2]/U[0] + U[7] * U[3]/U[0]);
		ret[5] = - U[1]/U[0] * U[6];
		ret[6] = 0.;
		ret[7] = - U[3]/U[0] * U[6];
	}
	else
		cout << __FUNCTION__ << ": err in component settings!";

	return ret;
}

double* NMgdMethods::psi_cusp(const double *U)//выделяет память под return массив!
{
	double *ret = new double[NInitial::getNcomp()];
	NMethods methods;

	ret[0] = 0.;
	ret[1] = 0.;
	ret[2] = 0.;
	ret[3] = 0.;
	ret[4] = methods.press(U) + .5 * (pow(U[5], 2) + pow(U[6], 2) + pow(U[7], 2));
	ret[5] = 0.;
	ret[6] = 0.;
	ret[7] = 0.;

	return ret;
}

void NMgdMethods::Roe(const double *Ur, const double *Ul, double *Flux, const char *dir, const char *dirn)
{
	double *f_r, *f_l, *P_l, *P_r, *psi_l, *psi_r;
	double M_l, M_r, a_half_l, a_half_r, a_half, C_plus, C_minus, alpha_plus_l, alpha_minus_r, beta_l;
	double beta_r, Dl_plus, Dr_minus, *psi_half;

	//double *traceData = new double[NInitial::getNcomp()];
	//for(int l = 0; l < NInitial::getNcomp(); l++)
	//	traceData[l] = Ur[l];
	//NTracer::traceToFile(traceData, "double", NInitial::getNcomp());
	//for(int l = 0; l < NInitial::getNcomp(); l++)
	//	traceData[l] = Ul[l];
	//NTracer::traceToFile(traceData, "double", NInitial::getNcomp());
	//delete []traceData;

	f_r = f_cusp(Ur, dir);
	f_l = f_cusp(Ul, dir);
	P_r = P_cusp(Ur, dir);
	P_l = P_cusp(Ul, dir);
	psi_r = psi_cusp(Ur);
	psi_l = psi_cusp(Ul);

	psi_half = new double[NInitial::getNcomp()];

	a_half_l = sound(Ul) + sound(Ul, dir);
	a_half_r = sound(Ur) + sound(Ur, dir);
	a_half = .5 * (a_half_l + a_half_r);
	if(!strcmp(dir, "r"))
	{
		M_l = Ul[1]/(Ul[0] * a_half);
		M_r = Ur[1]/(Ur[0] * a_half);
	}
	else
	{
		M_l = Ul[2]/(Ul[0] * a_half);
		M_r = Ur[2]/(Ur[0] * a_half);
	}
	alpha_plus_l = .5 * (1 + sign(M_l));
	alpha_minus_r = .5 * (1 - sign(M_r));
	beta_l = - max(0, 1 - static_cast< int >(abs(M_l)));
	beta_r = - max(0, 1 - static_cast< int >(abs(M_r)));
	C_plus = alpha_plus_l * (1 + beta_l) * M_l - .25 * beta_l * pow(M_l + 1, 2);
	C_minus = alpha_minus_r * (1 + beta_r) * M_r + .25 * beta_r * pow(M_r - 1, 2);
	Dl_plus = alpha_plus_l * (1 + beta_l) - .5 * beta_l * (1 + M_l);
	Dr_minus = alpha_minus_r * (1 + beta_r) - .5 * beta_r * (1 - M_r);
	
	for(size_t i = 0; i < NInitial::getNcomp(); i++)
	{
		psi_half[i] = a_half * (C_plus + C_minus) * (Dl_plus * psi_l[i] + Dr_minus * psi_r[i]);
		Flux[i] = a_half * (C_plus * f_l[i] + C_minus * f_r[i]) + Dl_plus * P_l[i] + Dr_minus * P_r[i] + 
		psi_half[i];
		if(Flux[i] != Flux[i])
		{
			cout << "Flux[" << i << "] is NAN!" << "\n";
			cout << "a_half_l " << a_half_l << "\n";
			cout << "a_half_r " << a_half_r << "\n";
			cout << "C_plus " << C_plus << "\n";
			cout << "C_minus " << C_minus << "\n";
			cout << "Dl_plus " << Dl_plus << "\n";
			cout << "Dr_minus " << Dr_minus << "\n";
			cout << "psi_l " << psi_l[i] << "\n";
			cout << "psi_r " << psi_r[i] << "\n";
			cout << "f_l " << f_l[i] << "\n";
			cout << "f_r " << f_r[i] << "\n";
			cout << "P_l " << P_l[i] << "\n";
			cout << "P_r " << P_r[i] << "\n";
			cout << "psi_half " << psi_half[i] << "\n";
			for(size_t i = 0; i < NInitial::getNcomp(); i++)
			{
				cout << "Ul[" << i << "] " << Ul[i] << "\n";
				cout << "Ur[" << i << "] " << Ur[i] << "\n";
			}
			cout << "sound(Ul) " << sound(Ul) << "\n";
			cout << "sound(Ur) " << sound(Ur) << "\n";
			cout << "sound(Ul, dir) " << sound(Ul, dir) << "\n";
			cout << "sound(Ur, dir) " << sound(Ur, dir) << "\n";
			cout << "press(Ur) " << press(Ur) << "\n";
			cout << "press(Ul) " << press(Ul) << "\n";

			throw 1;
		}
	}

	//double *traceData = new double[NInitial::getNcomp()];
	//for(int l = 0; l < NInitial::getNcomp(); l++)
	//	traceData[l] = Flux[l];
	//NTracer::traceToFile(traceData, "double", NInitial::getNcomp());
	//delete []traceData;

	delete []P_r;
	delete []P_l;

	delete []f_r;
	delete []f_l;

	delete []psi_r;
	delete []psi_l;

	delete []psi_half;

}
void NMgdMethods::Sourse(const double *U, double *S)
{
	NMethods::Sourse(U, S);
	S[1] += .5 * (pow(U[6], 2) + 3. * pow(U[7], 2) - pow(U[5], 2));
	S[2] -= U[5] * U[6];
	S[3] -= U[5] * U[7] + U[5] * U[6];
	S[4] += U[1]/U[0] * .5 * (pow(U[5], 2) + pow(U[6], 2) + pow(U[7], 2)) - 
		U[5] * (U[5] * U[1]/U[0] + U[6] * U[2]/U[0] + U[7] * U[3]/U[0]);
	S[5] = 0;
	S[6] = U[1]/U[0] * U[6] - U[2]/U[0] * U[5];
	S[7] = 2. * U[1]/U[0] * U[7];
}
vector<double> NMgdMethods::getSVariables(const double *Ur, const double *Ul)
{
	vector<double> sl, sr, s;
	sl.resize(8);
	sr.resize(8);
	s.resize(8);

	sl[0] = sqrt(Ul[0]);
	sl[1] = Ul[1]/sqrt(Ul[0]);
	sl[2] = Ul[2]/sqrt(Ul[0]);
	sl[3] = Ul[3]/sqrt(Ul[0]);
	sl[4] = (Ul[4] + .125/NInitial::get_pi() * (pow(Ul[5], 2) + pow(Ul[6], 2) + pow(Ul[7], 2)) + press(Ul))/sqrt(Ul[0]);
	sl[5] = Ul[5]/sqrt(Ul[0]);
	sl[6] = Ul[6]/sqrt(Ul[0]);
	sl[7] = Ul[7]/sqrt(Ul[0]);

	sr[0] = sqrt(Ur[0]);
	sr[1] = Ur[1]/sqrt(Ur[0]);
	sr[2] = Ur[2]/sqrt(Ur[0]);
	sr[3] = Ur[3]/sqrt(Ur[0]);
	sr[4] = (Ur[4] + .125/NInitial::get_pi() * (pow(Ur[5], 2) + pow(Ur[6], 2) + pow(Ur[7], 2)) + press(Ur))/sqrt(Ur[0]);
	sr[5] = Ur[5]/sqrt(Ur[0]);
	sr[6] = Ur[6]/sqrt(Ur[0]);
	sr[7] = Ur[7]/sqrt(Ur[0]);

	for(size_t i = 0; i < s.size(); i++)
		s[i] = .5 * (sr[i] + sl[i]);

	return s;
}
vector<vector<double> > NMgdMethods::roeMatrix(const double *Ur, const double *Ul, const char *dir, const char *dirn)
{
	double n[3];

	if(dir=="r")
	{
		if(dirn=="+")
			n[0]=1;
		else if(dirn=="-")
			n[0]=-1;
		else
		{
			cout<<"dirn in roe err";
			exit(1);
		}
		
		n[1]=0;
		n[2]=0;
	}
	else if(dir=="z")
	{
		if(dirn=="+")
			n[1]=1;
		else if(dirn=="-")
			n[1]=-1;
		else
		{
			cout<<"dirn in roe err";
			exit(1);
		}
		n[0]=0;
		n[2]=0;
	}
	else
	{
		cout<<"err in dir(Roe)"<<'\n';
		exit(1);
	}

	vector<vector<double> > A_Roe;
	vector<double> s(8);

	s = getSVariables(Ur, Ul);
	A_Roe.resize(8);
	for(size_t i = 0; i < 8; i++)
		A_Roe[i].resize(8);

	//filling the Roe matrix elements
	A_Roe[0][0] = 0.;
	A_Roe[0][1] = n[0];
	A_Roe[0][2] = n[1];
	A_Roe[0][3] = 0.;
	A_Roe[0][4] = 0.;
	A_Roe[0][5] = 0.;
	A_Roe[0][6] = 0.;
	A_Roe[0][7] = 0.;
	A_Roe[1][0] = -n[1] * s[2] * s[1]/pow(s[0], 2) + .5 * n[0] * (((pow(s[1], 2) + pow(s[2], 2) + pow(s[3], 2)) * NInitial::get_sigma() - 2 * pow(s[1], 2)) * 
		(NInitial::get_sigma() + 1) - NInitial::get_sigma() * 2 * pow(s[2], 2))/((1 + NInitial::get_sigma()) * pow(s[0], 2));
	A_Roe[1][1] = 1/s[0] * (s[2] * n[1] + s[1] * (2 - NInitial::get_sigma()) * n[0]);
	A_Roe[1][2] = 1/s[0] * (s[1] * n[1] - s[2] * NInitial::get_sigma() * n[0]);
	A_Roe[1][3] = -NInitial::get_sigma() * s[3]/s[0] * n[0];
	A_Roe[1][4] = NInitial::get_sigma() * n[0];
	A_Roe[1][5] = -s[0] * s[5] * ((NInitial::get_sigma() + 2)/(4 * NInitial::get_pi()) + 1) * n[0] - 1/(4 * NInitial::get_pi()) * s[0] * s[6] * n[1];
	A_Roe[1][6] = -.25 * s[0] * (s[5] * n[1] + s[6] * (4 + NInitial::get_sigma()/NInitial::get_pi()) * n[0]);
	A_Roe[1][7] = -s[0] * s[7] * (NInitial::get_sigma()/(4 * NInitial::get_pi()) + 1) * n[0];
	A_Roe[2][0] = n[1] * .5/pow(s[0], 2) * (NInitial::get_sigma() * (pow(s[1], 2) + pow(s[2], 2) + pow(s[3], 2)) - 2 * pow(s[2], 2)) -
		n[0] * s[1] * s[2]/pow(s[0], 2);
	A_Roe[2][1] = -1/s[0] * (s[1] * NInitial::get_sigma() * n[1] + s[2] * n[0]);
	A_Roe[2][2] = 1/s[0] * (n[1] * s[2] * (2 - NInitial::get_sigma()) + s[1] * n[0]);
	A_Roe[2][3] = -s[3] * NInitial::get_sigma()/s[0] * n[1];
	A_Roe[2][4] = NInitial::get_sigma() * n[1];
	A_Roe[2][5] = -n[1] * s[0] * s[5] * (NInitial::get_sigma() * (NInitial::get_sigma() - 3) * .25/((1 + NInitial::get_sigma()) * NInitial::get_pi()) + 1) - 
		s[0] * s[6] * n[0] * .25/NInitial::get_pi();
	A_Roe[2][6] = -s[0] * s[6] * n[1] * ((NInitial::get_sigma() + 2) * .25/NInitial::get_pi() - 1) - s[0] * s[5] * n[0] * .25/NInitial::get_pi();
	A_Roe[2][7] = n[1] * s[0] * s[7] * (1 - NInitial::get_sigma() * (NInitial::get_sigma() - 3) * .25/((1 + NInitial::get_sigma()) * NInitial::get_pi()));
	A_Roe[3][0] = -s[3]/pow(s[0], 2) * (s[2] * n[1] + s[1] * n[0]);
	A_Roe[3][1] = s[3]/s[0] * n[0];
	A_Roe[3][2] = s[3]/s[0] * n[1];
	A_Roe[3][3] = 1/s[0] * (s[2] * n[1] + s[1] * n[0]);
	A_Roe[3][4] = 0;
	A_Roe[3][5] = -.25 * s[0] * s[7] * n[0]/NInitial::get_pi();
	A_Roe[3][6] = -.25 * s[0] * s[7] * n[1]/NInitial::get_pi();
	A_Roe[3][7] = -.25 * s[0]/NInitial::get_pi() * (s[6] * n[1] + s[5] * n[0]);
	A_Roe[4][0] = n[1] * .5/s[0] * (s[2] * (-s[4] * (1 + NInitial::get_sigma())/s[0] + NInitial::get_sigma()/pow(s[0], 2) * (pow(s[1], 2) + pow(s[2], 2) + pow(s[3], 2))) + 
		.5 * s[6] * (s[5] * s[1] + s[7] * s[3] + s[6] * s[2])/NInitial::get_pi()) + n[0] * .5/pow(s[0], 2) * (s[1] * (pow(s[5], 2) + pow(s[6], 2) + pow(s[7], 2)) * (.25/NInitial::get_pi() + NInitial::get_sigma()/s[0]) -
		s[4] * s[1] * (1 + NInitial::get_sigma()) + .5 * pow(s[5] * s[0], 2)/NInitial::get_pi() * (-s[1] * (s[0] - 1)/pow(s[0], 2) * (1.5 * s[1] * s[5] + s[6] * s[2] + s[7] * s[3] - (pow(s[6], 2) - pow(s[7], 2))/pow(s[5], 2)) + 
		s[6] * s[2] + s[7] * s[3] + s[1] * s[5]));
	A_Roe[4][1] = -n[1] * (s[5] * s[6] * .25/NInitial::get_pi() + s[2] * s[1] * NInitial::get_sigma()/pow(s[0], 2)) + n[0] * (s[4]/s[0] + (pow(s[5], 2) + pow(s[6], 2) + pow(s[7], 2)) * (s[0] - 1) * .25 /(s[0] * NInitial::get_pi()) - 
		s[0] * pow(s[5], 3) * .25/NInitial::get_pi() - pow(s[1]/s[0], 2) * NInitial::get_sigma());
	A_Roe[4][2] = n[1] * (s[4]/s[0] - .25 * pow(s[6], 2)/NInitial::get_pi() - pow(s[2]/s[0], 2) * NInitial::get_sigma()) - n[0] * (.25 * s[0] * pow(s[5], 2) * s[6]/NInitial::get_pi() + s[1] * s[2] * NInitial::get_sigma()/pow(s[0], 2));
	A_Roe[4][3] = -n[1] * (.25 * s[7] * s[6]/NInitial::get_pi() + s[2] * s[3] * NInitial::get_sigma()/pow(s[0], 2)) - n[0] * (pow(s[5], 2) * s[0] * s[7] * .25/NInitial::get_pi() + s[1] * NInitial::get_sigma() * s[3]/pow(s[0], 2));
	A_Roe[4][4] = 1/s[0] * (n[1] * s[2] + s[1] * n[0]);
	A_Roe[4][5] = -.25 * n[1]/NInitial::get_pi() * (s[2] * s[5] * (NInitial::get_sigma() - 1) + s[6] * s[1]) + n[0] * s[5] * .5/NInitial::get_pi() * (-.5 * s[1] * (NInitial::get_sigma() - 1) + (s[1]/s[0] * (s[0] - 1) - s[0]) * 
		(s[6] * s[2] + s[7] * s[3] + 1.5 * s[5] * s[1]));
	A_Roe[4][6] = -n[1] * .25/NInitial::get_pi() * (s[2] * s[6] * (NInitial::get_sigma() + 1) + s[5] * s[1] + s[7] * s[3]) + n[0] * .5/NInitial::get_pi() * (-.5 * (s[1] * s[6] * (NInitial::get_sigma() - 1) + s[0] * pow(s[5], 2) * s[2]) + 
		1/s[0] * (s[0] - 1) * s[6] * s[1]);
	A_Roe[4][7] = -.25 * n[1]/NInitial::get_pi() * (s[2] * s[7] * (NInitial::get_sigma() - 1) + s[6] * s[3]) + n[0] * .25/(s[0] * NInitial::get_pi()) * (s[1] * s[7] * (s[0] * (NInitial::get_sigma() + 1) - 2) - pow(s[5] * s[0], 2) * s[3]);
	A_Roe[5][0] = n[1]/pow(s[0], 2) * (s[1] * s[6] - s[5] * s[2]);
	A_Roe[5][1] = -n[1] * s[6]/s[0];
	A_Roe[5][2] = n[1] * s[5]/s[0];
	A_Roe[5][3] = 0.;
	A_Roe[5][4] = 0.;
	A_Roe[5][5] = s[2] * n[1]/s[0];
	A_Roe[5][6] = -s[1] * n[1]/s[0];
	A_Roe[5][7] = 0.;
	A_Roe[6][0] = n[0]/pow(s[0], 2) * (s[5] * s[2] - s[6] * s[1]);
	A_Roe[6][1] = s[6]/s[0] * n[0];
	A_Roe[6][2] = -s[5]/s[0] * n[0];
	A_Roe[6][3] = 0.;
	A_Roe[6][4] = 0.;
	A_Roe[6][5] = -s[2] * n[0]/s[0];
	A_Roe[6][6] = s[1] * n[0]/s[0];
	A_Roe[6][7] = 0.;
	A_Roe[7][0] = n[1]/pow(s[0], 2) * (s[6] * s[3] - s[7] * s[2]) + n[0] * s[5] * s[3]/pow(s[0], 2);
	A_Roe[7][1] = s[7]/s[0] * n[0];
	A_Roe[7][2] = s[7]/s[0] * n[1];
	A_Roe[7][3] = -(s[6] * n[1] + s[5] * n[0])/s[0];
	A_Roe[7][4] = 0.;
	A_Roe[7][5] = -s[3] * n[0]/s[0];
	A_Roe[7][6] = -s[3] * n[1]/s[0];
	A_Roe[7][7] = (s[2] * n[1] + s[1] * n[0])/s[0];

	return A_Roe;
}
double NMgdMethods::energy(const double *U)
{
	return (U[4] - .5 * (pow(U[5], 2) + pow(U[6], 2) + pow(U[7], 2)))/U[0] 
		- .5 * (pow(U[1]/U[0], 2) + pow(U[2]/U[0], 2) + pow(U[3]/U[0], 2));
}

double*** NMgdMethods::korrectH(const double ***U, const double *hr, const double *hz, double &tau)//выделяет память под return массив!
{
	double ***H_Korrect;//[r,z,phi]
	double *F;
	double Fz_i_jmh, Fz_im_jmh, Fr_imh_jm, Fz_i_jph, Fz_im_jph, Fr_imh_jp, Fz_ip_jmh, Fr_iph_j, Fr_iph_jm, Fr_imh_j;
	NMgdBoundaryCond bound;
	vector<int> sizes(2);
	sizes[0] = NInitial::get_xmax();
	sizes[1] = NInitial::get_ymax();
	NArrPacker<double> *U_pack = new NArrPacker<double>(DUMMY_NUM, sizes, NInitial::getPeriodicName(), NInitial::getSlipName(), (void*)U);
	double*** U_pack_arr = (double***)U_pack->getPackArr();

	F = new double[NInitial::getNcomp()];

	H_Korrect = new double**[NInitial::get_xmax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		H_Korrect[i] = new double*[NInitial::get_ymax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			H_Korrect[i][j] = new double[3];

	//double *traceData = new double[10];
	//for(int i = 0; i <= NInitial::get_xmax() + 1; i++)
	//	for(int j = 0; j <= NInitial::get_ymax() + 1; j++)
	//	{
	//		traceData[0] = U_pack_arr[i][j][0];
	//		traceData[1] = U_pack_arr[i][j][1];
	//		traceData[2] = U_pack_arr[i][j][2];
	//		traceData[3] = U_pack_arr[i][j][3];
	//		traceData[4] = U_pack_arr[i][j][4];
	//		traceData[5] = U_pack_arr[i][j][5];
	//		traceData[6] = U_pack_arr[i][j][6];
	//		traceData[7] = U_pack_arr[i][j][7];

	//		traceData[8] = i;
	//		traceData[9] = j;

	//		NTracer::traceToFile(traceData, "double", 10);
	//	}

	//	exit(1);

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
		{
			Roe(U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i + DUMMY_NUM][j - 1 + DUMMY_NUM], F, "z", NULL);
			Fz_i_jmh = F[5];
			Roe(U_pack_arr[i - 1 + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i - 1 + DUMMY_NUM][j - 1 + DUMMY_NUM], F, "z", NULL);
			Fz_im_jmh = F[5];
			Roe(U_pack_arr[i + DUMMY_NUM][j - 1 + DUMMY_NUM], U_pack_arr[i - 1 + DUMMY_NUM][j - 1 + DUMMY_NUM], F, "r", NULL);
			Fr_imh_jm = F[6];
			Roe(U_pack_arr[i + DUMMY_NUM][j + 1 + DUMMY_NUM], U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], F, "z", NULL);
			Fz_i_jph = F[5];
			Roe(U_pack_arr[i - 1 + DUMMY_NUM][j + 1 + DUMMY_NUM], U_pack_arr[i - 1 + DUMMY_NUM][j + DUMMY_NUM], F, "z", NULL);
			Fz_im_jph = F[5];
			Roe(U_pack_arr[i + DUMMY_NUM][j + 1 + DUMMY_NUM], U_pack_arr[i - 1 + DUMMY_NUM][j + 1 + DUMMY_NUM], F, "r", NULL);
			Fr_imh_jp = F[6];
			Roe(U_pack_arr[i + 1 + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i + 1 + DUMMY_NUM][j - 1 + DUMMY_NUM], F, "z", NULL);
			Fz_ip_jmh = F[5];
			Roe(U_pack_arr[i + 1 + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], F, "r", NULL);
			Fr_iph_j = F[6];
			Roe(U_pack_arr[i + 1 + DUMMY_NUM][j - 1 + DUMMY_NUM], U_pack_arr[i + DUMMY_NUM][j - 1 + DUMMY_NUM], F, "r", NULL);
			Fr_iph_jm = F[6];
			Roe(U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i - 1 + DUMMY_NUM][j + DUMMY_NUM], F, "r", NULL);
			Fr_imh_j = F[6];

			H_Korrect[i][j][0] = U[i][j][5] + tau * .25/hz[j] * (Fz_i_jmh + Fz_im_jmh - Fr_imh_jm - Fz_i_jph 
				- Fz_im_jph + Fr_imh_jp);
			H_Korrect[i][j][1] = U[i][j][6] + tau * .25/hr[i] * (Fz_ip_jmh - Fr_iph_j - Fr_iph_jm - Fz_im_jmh + 
				Fr_imh_j + Fr_imh_jm);


			Roe(U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i + DUMMY_NUM][j - 1 + DUMMY_NUM], F, "z", NULL);
			Fz_i_jmh = F[7];
			Roe(U_pack_arr[i + DUMMY_NUM][j + 1 + DUMMY_NUM], U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], F, "z", NULL);
			Fz_i_jph = F[7];
			Roe(U_pack_arr[i + 1 + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], F, "r", NULL);
			Fr_iph_j = F[7];
			Roe(U_pack_arr[i + DUMMY_NUM][j + DUMMY_NUM], U_pack_arr[i - 1 + DUMMY_NUM][j + DUMMY_NUM], F, "r", NULL);
			Fr_imh_j = F[7];

			H_Korrect[i][j][2] = U[i][j][7] + tau * (1/hr[i] * (Fr_imh_j - Fr_iph_j) + 1/hz[j] * (Fz_i_jmh - Fz_i_jph));

			//traceData[0] = H_Korrect[i][j][0];
			//traceData[1] = H_Korrect[i][j][1];
			//traceData[2] = H_Korrect[i][j][2];

			//traceData[3] = i;
			//traceData[4] = j;

			//NTracer::traceToFile(traceData, "double", 5);
		}


	//delete []traceData;

	delete []F;
	delete U_pack;

	return H_Korrect;
}

