#pragma once
#ifndef METHODS_CPP
#define METHODS_CPP

#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\methods.h"
#else
#include "includes/methods.h"
#endif

NMethods::~NMethods()
{
	;
}
//vector<double> NMethods::getSVariables(const double *Ur, const double *Ul)
//{
//	vector<double> sl, sr, s;
//	sl.resize(NInitial::getNcomp());
//	sr.resize(NInitial::getNcomp());
//	s.resize(NInitial::getNcomp());
//
//	sl[0] = sqrt(Ul[0]);
//	sl[1] = Ul[1]/sqrt(Ul[0]);
//	sl[2] = Ul[2]/sqrt(Ul[0]);
//	sl[3] = Ul[3]/sqrt(Ul[0]);
//	sl[4] = (Ul[4] + press(Ul))/sqrt(Ul[0]);
//
//	sr[0] = sqrt(Ur[0]);
//	sr[1] = Ur[1]/sqrt(Ur[0]);
//	sr[2] = Ur[2]/sqrt(Ur[0]);
//	sr[3] = Ur[3]/sqrt(Ur[0]);
//	sr[4] = (Ur[4] + press(Ur))/sqrt(Ur[0]);
//
//	for(size_t i = 0; i < s.size(); i++)
//		s[i] = .5 * (sr[i] + sl[i]);
//
//	return s;
//}
//vector<vector<double> > NMethods::roeMatrix(const double *Ur, const double *Ul, const char *dir, const char *dirn)
//{
//	double n[3];
//
//	if(dir=="r")
//	{
//		if(dirn=="+")
//			n[0]=1;
//		else if(dirn=="-")
//			n[0]=-1;
//		else
//		{
//			cout<<"dirn in roe err";
//			exit(1);
//		}
//		
//		n[1]=0;
//		n[2]=0;
//	}
//	else if(dir=="z")
//	{
//		if(dirn=="+")
//			n[1]=1;
//		else if(dirn=="-")
//			n[1]=-1;
//		else
//		{
//			cout<<"dirn in roe err";
//			exit(1);
//		}
//		n[0]=0;
//		n[2]=0;
//	}
//	else
//	{
//		cout<<"err in dir(Roe)"<<'\n';
//		exit(1);
//	}
//
//	vector<vector<double> > A_Roe;
//	vector<double> s(NInitial::getNcomp());
//
//	s = getSVariables(Ur, Ul);
//	A_Roe.resize(NInitial::getNcomp());
//	for(size_t i = 0; i < NInitial::getNcomp(); i++)
//		A_Roe[i].resize(NInitial::getNcomp());
//
//	//filling the Roe matrix elements
//	A_Roe[0][0] = 0.;
//	A_Roe[0][1] = n[0];
//	A_Roe[0][2] = n[1];
//	A_Roe[0][3] = 0.;
//	A_Roe[0][4] = 0.;
//	A_Roe[1][0] = -n[1] * s[2] * s[1]/pow(s[0], 2) + .5 * n[0] * (((pow(s[1], 2) + pow(s[2], 2) + pow(s[3], 2)) * NInitial::get_sigma() - 2 * pow(s[1], 2)) * 
//		(NInitial::get_sigma() + 1) - NInitial::get_sigma() * 2 * pow(s[2], 2))/((1 + NInitial::get_sigma()) * pow(s[0], 2));
//	A_Roe[1][1] = 1/s[0] * (s[2] * n[1] + s[1] * (2 - NInitial::get_sigma()) * n[0]);
//	A_Roe[1][2] = 1/s[0] * (s[1] * n[1] - s[2] * NInitial::get_sigma() * n[0]);
//	A_Roe[1][3] = -NInitial::get_sigma() * s[3]/s[0] * n[0];
//	A_Roe[1][4] = NInitial::get_sigma() * n[0];
//	A_Roe[2][0] = n[1] * .5/pow(s[0], 2) * (NInitial::get_sigma() * (pow(s[1], 2) + pow(s[2], 2) + pow(s[3], 2)) - 2 * pow(s[2], 2)) -
//		n[0] * s[1] * s[2]/pow(s[0], 2);
//	A_Roe[2][1] = -1/s[0] * (s[1] * NInitial::get_sigma() * n[1] + s[2] * n[0]);
//	A_Roe[2][2] = 1/s[0] * (n[1] * s[2] * (2 - NInitial::get_sigma()) + s[1] * n[0]);
//	A_Roe[2][3] = -s[3] * NInitial::get_sigma()/s[0] * n[1];
//	A_Roe[2][4] = NInitial::get_sigma() * n[1];
//	A_Roe[3][0] = -s[3]/pow(s[0], 2) * (s[2] * n[1] + s[1] * n[0]);
//	A_Roe[3][1] = s[3]/s[0] * n[0];
//	A_Roe[3][2] = s[3]/s[0] * n[1];
//	A_Roe[3][3] = 1/s[0] * (s[2] * n[1] + s[1] * n[0]);
//	A_Roe[3][4] = 0;
//	A_Roe[4][0] = n[1] * .5/s[0] * (s[2] * (-s[4] * (1 + NInitial::get_sigma())/s[0] + NInitial::get_sigma()/pow(s[0], 2) * (pow(s[1], 2) + pow(s[2], 2) + pow(s[3], 2)))) + n[0] * .5/pow(s[0], 2) * (-s[4] * s[1] * (1 + NInitial::get_sigma()));
//	A_Roe[4][1] = -n[1] * (s[2] * s[1] * NInitial::get_sigma()/pow(s[0], 2)) + n[0] * (s[4]/s[0] + (pow(s[5], 2) + pow(s[6], 2) + pow(s[7], 2)) * (s[0] - 1) * .25 /(s[0] * NInitial::get_pi()) - 
//		s[0] * pow(s[5], 3) * .25/NInitial::get_pi() - pow(s[1]/s[0], 2) * NInitial::get_sigma());
//	A_Roe[4][2] = n[1] * (s[4]/s[0] - .25 * pow(s[6], 2)/NInitial::get_pi() - pow(s[2]/s[0], 2) * NInitial::get_sigma()) - n[0] * (.25 * s[0] * pow(s[5], 2) * s[6]/NInitial::get_pi() + s[1] * s[2] * NInitial::get_sigma()/pow(s[0], 2));
//	A_Roe[4][3] = -n[1] * (.25 * s[7] * s[6]/NInitial::get_pi() + s[2] * s[3] * NInitial::get_sigma()/pow(s[0], 2)) - n[0] * (pow(s[5], 2) * s[0] * s[7] * .25/NInitial::get_pi() + s[1] * NInitial::get_sigma() * s[3]/pow(s[0], 2));
//	A_Roe[4][4] = 1/s[0] * (n[1] * s[2] + s[1] * n[0]);
//	A_Roe[4][5] = -.25 * n[1]/NInitial::get_pi() * (s[2] * s[5] * (NInitial::get_sigma() - 1) + s[6] * s[1]) + n[0] * s[5] * .5/NInitial::get_pi() * (-.5 * s[1] * (NInitial::get_sigma() - 1) + (s[1]/s[0] * (s[0] - 1) - s[0]) * 
//		(s[6] * s[2] + s[7] * s[3] + 1.5 * s[5] * s[1]));
//	A_Roe[4][6] = -n[1] * .25/NInitial::get_pi() * (s[2] * s[6] * (NInitial::get_sigma() + 1) + s[5] * s[1] + s[7] * s[3]) + n[0] * .5/NInitial::get_pi() * (-.5 * (s[1] * s[6] * (NInitial::get_sigma() - 1) + s[0] * pow(s[5], 2) * s[2]) + 
//		1/s[0] * (s[0] - 1) * s[6] * s[1]);
//	A_Roe[4][7] = -.25 * n[1]/NInitial::get_pi() * (s[2] * s[7] * (NInitial::get_sigma() - 1) + s[6] * s[3]) + n[0] * .25/(s[0] * NInitial::get_pi()) * (s[1] * s[7] * (s[0] * (NInitial::get_sigma() + 1) - 2) - pow(s[5] * s[0], 2) * s[3]);
//	A_Roe[5][0] = n[1]/pow(s[0], 2) * (s[1] * s[6] - s[5] * s[2]);
//	A_Roe[5][1] = -n[1] * s[6]/s[0];
//	A_Roe[5][2] = n[1] * s[5]/s[0];
//	A_Roe[5][3] = 0.;
//	A_Roe[5][4] = 0.;
//	A_Roe[5][5] = s[2] * n[1]/s[0];
//	A_Roe[5][6] = -s[1] * n[1]/s[0];
//	A_Roe[5][7] = 0.;
//	A_Roe[6][0] = n[0]/pow(s[0], 2) * (s[5] * s[2] - s[6] * s[1]);
//	A_Roe[6][1] = s[6]/s[0] * n[0];
//	A_Roe[6][2] = -s[5]/s[0] * n[0];
//	A_Roe[6][3] = 0.;
//	A_Roe[6][4] = 0.;
//	A_Roe[6][5] = -s[2] * n[0]/s[0];
//	A_Roe[6][6] = s[1] * n[0]/s[0];
//	A_Roe[6][7] = 0.;
//	A_Roe[7][0] = n[1]/pow(s[0], 2) * (s[6] * s[3] - s[7] * s[2]) + n[0] * s[5] * s[3]/pow(s[0], 2);
//	A_Roe[7][1] = s[7]/s[0] * n[0];
//	A_Roe[7][2] = s[7]/s[0] * n[1];
//	A_Roe[7][3] = -(s[6] * n[1] + s[5] * n[0])/s[0];
//	A_Roe[7][4] = 0.;
//	A_Roe[7][5] = -s[3] * n[0]/s[0];
//	A_Roe[7][6] = -s[3] * n[1]/s[0];
//	A_Roe[7][7] = (s[2] * n[1] + s[1] * n[0])/s[0];
//
//	return A_Roe;
//}
//
double NMethods::sign(const double &val)
{
	return val > 0 ? 1 : (val < 0 ? -1 : 0);
}
void NMethods::LoadGrid(double *hr, double *hz)
{//Si = h[0] * (1-pow(RFAKTOR,i))/(1-RFAKTOR), i=[0,n]
	//hr[0] = NInitial::get_deltar()/(2 * (1 - pow(NInitial::getRFAKTOR(), NInitial::getRZONE() - 1))/(1 - NInitial::getRFAKTOR()) + pow(NInitial::getRFAKTOR(), NInitial::getRZONE() - 1) * (NInitial::get_xmax() - 2 * 
	//	NInitial::getRZONE()));
	hr[0]=NInitial::get_deltar()/(NInitial::get_xmax());
	for(int i = 0; i < NInitial::get_xmax() - 1; i++)
	{
		if((i + 1) < NInitial::getRZONE())
			hr[i + 1] = hr[i] * NInitial::getRFAKTOR();
		else if((i + 1) > (NInitial::get_xmax() - NInitial::getRZONE()))
			hr[i + 1]=hr[i]/NInitial::getRFAKTOR();
		else
			hr[i + 1]=hr[i];
	}
	
	hr[NInitial::get_xmax()] = hr[NInitial::get_xmax() - 1];
	hr[NInitial::get_xmax() + 1] = hr[NInitial::get_xmax() - 2];
	hr[NInitial::get_xmax() + 2] = hr[NInitial::get_xmax() - 3];

	for(int j = 0; j < NInitial::get_ymax() + 2; j++)
		hz[j] = NInitial::get_deltaz()/NInitial::get_ymax();

	hz[NInitial::get_ymax()] = hz[NInitial::get_ymax() - 1];
	hz[NInitial::get_ymax() + 1] = hz[NInitial::get_ymax() - 2];
	hz[NInitial::get_ymax() + 2] = hz[NInitial::get_ymax() - 3];
}
double NMethods::zij1(const int &j, const double &h0)
{
	if(j >= 0)
		return h0 * j + h0 * .5;
	else
		return h0 * (abs(j) - 1) + h0 * .5;
}
double NMethods::press(const double *U)
{	
	double V;
	V = pow(U[1]/U[0], 2) + pow(U[2]/U[0], 2) + pow(U[3]/U[0], 2);
	if(U[3] != U[3])
	{
		cout << __FUNCTION__ << " U[4]: NAN, " << "U[0]=" << U[0] << ",U[1]=" << U[1] << ",U[2]=" << U[2]
		 << ",U[3]=" << U[3] << ",U[4]=" << U[4] << ",U[5]=" << U[5] << ",U[6]=" << U[6];
		exit(1);
	}
	else if(V != V)
	{
		cout << __FUNCTION__ << "V: NAN";
		exit(1);
	}

	return NInitial::get_sigma() * U[0] * (U[4]/U[0] - V * .5 - .5/U[0] * (pow(U[5], 2) + pow(U[6], 2) + pow(U[7], 2))); //additional term in the MGD case
}
double NMethods::sound(const double *U)
{
	double p, ro;
	const double eps = .000001;

	p = press(U);
	ro = U[0];
	if(p != p | ro != ro)
	{
		cout << __FUNCTION__ << "p or ro: NAN or INF!";
		for(int k=0; k < NInitial::getNcomp(); k++)
		{
			cout << "U[" << k << "] = " << U[k] << "\n";
		}
		exit(1);
	}

	return sqrt((NInitial::get_sigma() + 1) * p/ro);
}
double NMethods::MakhNumber(const double *U)
{
	return sqrt(pow(U[1]/U[0], 2) + pow(U[2]/U[0], 2)/* + pow(U[3]/U[0], 2)*/)/sound(U);
}
double NMethods::entalpy(const double *U)
{
	return U[3]/U[0] + press(U)/U[0];
}
void NMethods::get_charact_variables(const double *U, double *V)
{
	V[0] = U[0];
	V[1] = U[1]/U[0];
	V[2] = U[2]/U[0];
	//V[3] = U[3]/U[0];
	V[3] = press(U);
}
void NMethods::Frp(const double *U, double *F)
{
	F[0] = U[1];
	F[1] = pow(U[1], 2)/U[0] + press(U);
	F[2] = U[2] * U[1]/U[0];
	F[3] = U[1] * U[3]/U[0];
	F[4] = (U[3] + press(U)) * U[1]/U[0];
}
void NMethods::Fz(const double *U, double *F)
{
	F[0] = U[2];
	F[1] = U[2] * U[1]/U[0];
	F[2] = pow(U[2], 2)/U[0] + press(U);
	F[3] = U[2] * U[3]/U[0];
	F[4] = (U[3] + press(U)) * U[2]/U[0];
}
void NMethods::Roe(const double *Ur, const double *Ul, double *Flux, const char *dir, const char *dirn)
{
	const double deltaentropy = 0.374165739;
	double *Vr, *Vl, El, Er, Hr, Hl, rotilda, vtilda[3], Htilda, qtilda, ctilda, n[3], Vtilda,
		vl[3], vr[3], deltap, deltaro, deltaV, deltau, deltav, deltaw, Vect1[5], Vect2[5], Vect3[5],
		Vect4[5], deltaF1[5], deltaF234[5], deltaF5[5], *Fr, *Fl, Velr, Vell/*,n1[4]*/;

	Vr = new double[NInitial::getNcomp()];
	Vl = new double[NInitial::getNcomp()];
	Fr = new double[NInitial::getNcomp()];
	Fl = new double[NInitial::getNcomp()];

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
		Frp(Ur,Fr);
		Frp(Ul,Fl);
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
		Fz(Ur,Fr);
		Fz(Ul,Fl);
	}
	else
	{
		cout<<"err in dir(Roe)"<<'\n';
		exit(1);
	}
		
	get_charact_variables(Ur,Vr);
	get_charact_variables(Ul,Vl);

	El=Ul[4]/Ul[0];
	Er=Ur[4]/Ur[0];
	Hr=entalpy(Ur);
	Hl=entalpy(Ul);

	rotilda=sqrt(Vr[0]*Vl[0]);
	vtilda[0]=(Vl[1]*sqrt(Vl[0])+Vr[1]*sqrt(Vr[0]))/(sqrt(Vl[0])+sqrt(Vr[0]));
	vtilda[1]=(Vl[2]*sqrt(Vl[0])+Vr[2]*sqrt(Vr[0]))/(sqrt(Vl[0])+sqrt(Vr[0]));
	vtilda[2]=(Vl[3]*sqrt(Vl[0])+Vr[3]*sqrt(Vr[0]))/(sqrt(Vl[0])+sqrt(Vr[0]));
	Htilda=(Hl*sqrt(Vl[0])+Hr*sqrt(Vr[0]))/(sqrt(Vl[0])+sqrt(Vr[0]));
	qtilda=pow(vtilda[0],2)+pow(vtilda[1],2)+pow(vtilda[2],2);
	ctilda=sqrt(NInitial::get_sigma()*(Htilda-qtilda*.5));

		//get_VandS_inf(N,Nsurf,V,S,n,a,hr,hz,hphi,axiscoord);

	Vtilda=vtilda[0]*n[0]+vtilda[1]*n[1]+vtilda[2]*n[2];
	vl[0]=Vl[1];
	vl[1]=Vl[2];
	vl[2]=Vl[3];
	vr[0]=Vr[1];
	vr[1]=Vr[2];
	vr[2]=Vr[3];
	Velr=vr[0]*n[0]+vr[1]*n[1]+vr[2]*n[2];//a.get_scalarmultiply(vr,n);//n1
	Vell=vl[0]*n[0]+vl[1]*n[1]+vl[2]*n[2];//a.get_scalarmultiply(vl,n);//n1

	deltap=Vr[4]-Vl[4];
	deltaV=Velr-Vell;
	deltaro=Vr[0]-Vl[0];
	deltau=vr[0]-vl[0];
	deltav=vr[1]-vl[1];
	deltaw=vr[2]-vl[2];

	Vect1[0]=1;
	Vect1[1]=vtilda[0]-ctilda*n[0];//what is normal (outer or insider)?
	Vect1[2]=vtilda[1]-ctilda*n[1];//what is normal (outer or insider)?
	Vect1[3]=vtilda[2]-ctilda*n[2];//what is normal (outer or insider)?
	Vect1[4]=Htilda-ctilda*Vtilda;

	Vect2[0]=1;
	Vect2[1]=vtilda[0];
	Vect2[2]=vtilda[1];
	Vect2[3]=vtilda[2];
	Vect2[4]=qtilda*.5;
		
	Vect3[0]=0;
	Vect3[1]=deltau-deltaV*n[0];//what is normal (outer or insider)?
	Vect3[2]=deltav-deltaV*n[1];//what is normal (outer or insider)?
	Vect3[3]=deltaw-deltaV*n[2];//what is normal (outer or insider)?
	Vect3[4]=vtilda[0]*deltau+vtilda[1]*deltav+vtilda[2]*deltaw-deltaV*Vtilda;

	Vect4[0]=1;
	Vect4[1]=vtilda[0]+ctilda*n[0];//what is normal (outer or insider)?
	Vect4[2]=vtilda[1]+ctilda*n[1];//what is normal (outer or insider)?
	Vect4[3]=vtilda[2]+ctilda*n[2];//what is normal (outer or insider)?
	Vect4[4]=Htilda+ctilda*Vtilda;

	for(int i=0;i<=4;i++)
	{
		deltaF1[i]=.5*abs(Vtilda-ctilda)*((deltap-rotilda*ctilda*deltaV)/pow(ctilda,2))*Vect1[i];
		deltaF234[i]=abs(Vtilda)*((deltaro-deltap/pow(ctilda,2))*Vect2[i]+rotilda*Vect3[i]);
		deltaF5[i]=.5*abs(Vtilda+ctilda)*((deltap+rotilda*ctilda*deltaV)/pow(ctilda,2))*Vect4[i];

		if(abs(Vtilda-ctilda) <= deltaentropy)
	     	deltaF1[i]=.5*((pow(Vtilda-ctilda,2)+pow(deltaentropy,2))*.5/deltaentropy)*((deltap-rotilda*ctilda*deltaV)/pow(ctilda,2))*Vect1[i];

		if(abs(Vtilda+ctilda)<=deltaentropy)
			deltaF5[i]=.5*((pow(Vtilda+ctilda,2)+pow(deltaentropy,2))*.5/deltaentropy)*((deltap-rotilda*ctilda*deltaV)/pow(ctilda,2))*Vect1[i];
//if(((deltaF1[i]!=0) || (deltaF234[i]!=0) || (deltaF5[i]!=0)))// && (B[Nsurf]!=1))
//cout<<"a"<<' '<<N<<' ';
		Flux[i]=.5*(Fr[i]+Fl[i]-deltaF1[i]-deltaF234[i]-deltaF5[i]);
//if(((Fr[i]+Fl[i])!=0) && (B[Nsurf]!=1))
//cout<<"a1"<<' '<<N<<' ';

	}
	
	delete []Vr;
	delete []Vl;
	delete []Fr;
	delete []Fl;
}
void NMethods::Sourse(const double *U, double *S)
{
	S[0]=U[1];
	S[1]=pow(U[1],2)/U[0]-pow(U[3],2)/U[0];
	S[2]=U[2]*U[1]/U[0];
	S[3]=2*U[1]*U[3]/U[0];
	S[4]=(U[4]+press(U))*U[1]/U[0];
	//S[0]=U[1];
	//S[1]=(pow(U[1],2)-pow(U[3],2))/U[0];
	//S[2]=U[1]*U[2]/U[0];
	//S[3]=2*U[1]*U[3]/U[0];
	//S[4]=(U[4]+press(U))*U[1]/U[0];
}
void NMethods::get_grad1(const int i, const int j, double ***U, double **grad, const char *bz, const char *br,
	const double *hr, const double *hz)
{
	double Sr1[2], Sr2[2], Sz[2], *Uzp1, *Urp1, lambdafrp, lambdafzm, lambdafrm, lambdafzp, vol,
		*Uzm1, *Urm1;

	Uzp1 = new double[NInitial::getNcomp()];
	Urp1 = new double[NInitial::getNcomp()];
	Uzm1 = new double[NInitial::getNcomp()];
	Urm1 = new double[NInitial::getNcomp()];

	if((i>=0) && (j>=0))
	{
		Sr1[0]=hz[j];
		Sr2[0]=hz[j];
		Sz[1]=hr[i];
		vol=hz[j]*hr[i];
	}
	else if((i<0) && (j<0))
	{
		Sr1[0]=hz[abs(j)-1];
		Sr2[0]=hz[abs(j)-1];
		Sz[1]=hr[abs(i)-1];
		vol=hz[abs(j)-1]*hr[abs(i)-1];
	}
	else if((i<0) && (j>=0))
	{
		Sr1[0]=hz[j];
		Sr2[0]=hz[j];
		Sz[1]=hr[abs(i)-1];
		vol=hz[j]*hr[abs(i)-1];
	}
	else if((j<0) && (i>=0))
	{
		Sz[1]=hr[i];
		Sr1[0]=hz[abs(j)-1];
		Sr2[0]=hz[abs(j)-1];
		vol=hz[abs(j)-1]*hr[i];
	}

	Sr1[1]=0;
	Sr2[1]=0;
	Sz[0]=0;

	if(i>0)
	{
		lambdafrp=hr[i]/(hr[i]+hr[i+1]);
		lambdafrm=hr[i]/(hr[i]+hr[i-1]);
	}
	else if(i==0)
	{
		lambdafrp=hr[i]/(hr[i]+hr[i+1]);
		lambdafrm=hr[i]/(hr[i]+hr[0]);
	}
	else if(i==-1)
	{
		lambdafrp=hr[abs(i)-1]/(hr[abs(i)-1]+hr[i+1]);
		lambdafrm=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i-1)-1]);
	}
	else 
	{
		lambdafrp=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i+1)-1]);
		lambdafrm=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i-1)-1]);
	}

	if(j>0)
	{
		lambdafzp=hz[j]/(hz[j]+hz[j+1]);
		lambdafzm=hz[j]/(hz[j]+hz[j-1]);
	}
	else if(j==0)
	{
		lambdafzp=hz[j]/(hz[j]+hz[j+1]);
		lambdafzm=hz[j]/(hz[j]+hz[abs(j-1)-1]);
	}
	else if(j==-1)
	{
		lambdafzp=hz[abs(j)-1]/(hz[abs(j)-1]+hz[j+1]);
		lambdafzm=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j-1)-1]);
	}
	else 
	{
		lambdafzp=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j+1)-1]);
		lambdafzm=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j-1)-1]);
	}

	double wxp1,wxm1,wyp1,wym1,deltarxp1,deltarym1,deltarxm1,deltaryp1;

	if(i>-1)
	{
		if(j>0)
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[j]+hz[j-1])*.5);
		}
		else if(j==0)
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[j]+hz[abs(j-1)-1])*.5);
		}
		else if(j==-1)
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[abs(j)-1]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
		else
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[abs(j)-1]+hz[abs(j+1)-1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
	}
	else
	{
		if(j>0)
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[j]+hz[j-1])*.5);
		}
		else if(j==0)
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[j]+hz[abs(j-1)-1])*.5);
		}
		else if(j==-1)
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
		else
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[abs(j+1)-1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
	}


	if(j>-1)
	{
		if(i>0)
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[i]+hr[i-1])*.5);
		}
		else if(i==0)
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[i]+hr[abs(i-1)-1])*.5);
		}
		else if(i==-1)
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[abs(i)-1]+hr[i+1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
		else
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[abs(i)-1]+hr[abs(i+1)-1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
	}
	else
	{
		if(i>0)
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[i]+hr[i-1])*.5);
		}
		else if(i==0)
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[i]+hr[abs(i-1)-1])*.5);
		}
		else if(i==-1)
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[i+1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
		else
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[abs(i+1)-1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
	}
				
	if(j>0)
	{
		deltarxp1=(hz[j]+hz[j+1])*.5;
		deltarxm1=(hz[j]+hz[j-1])*.5;
	}
	else if(j==0)
	{
		deltarxp1=(hz[j]+hz[j+1])*.5;
		deltarxm1=(hz[j]+hz[abs(j-1)-1])*.5;
	}
	else if(j==-1)
	{
		deltarxp1=(hz[abs(j)-1]+hz[j+1])*.5;
		deltarxm1=(hz[abs(j)-1]+hz[abs(j-1)-1])*.5;
	}
	else
	{
		deltarxp1=(hz[abs(j)-1]+hz[abs(j+1)-1])*.5;
		deltarxm1=(hz[abs(j)-1]+hz[abs(j-1)-1])*.5;
	}

	if(i>0)
	{
		deltaryp1=(hr[i]+hr[i+1])*.5;
		deltarym1=(hr[i]+hr[i-1])*.5;
	}
	else if(i==0)
	{
		deltaryp1=(hr[i]+hr[i+1])*.5;
		deltarym1=(hr[i]+hr[abs(i-1)-1])*.5;
	}
	else if(i==-1)
	{
		deltaryp1=(hr[abs(i)-1]+hr[i+1])*.5;
		deltarym1=(hr[abs(i)-1]+hr[abs(i-1)-1])*.5;
	}
	else
	{
		deltaryp1=(hr[abs(i)-1]+hr[abs(i+1)-1])*.5;
		deltarym1=(hr[abs(i)-1]+hr[abs(i-1)-1])*.5;
	}
	
	for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
	{
		Uzp1[k]=/*U[i][j][k]*(1-lambdafzp)+*/U[i + DUMMY_NUM2][j + 1 + DUMMY_NUM2][k]/**lambdafzp*/;

		Uzm1[k]=/*U[i][j][k]*(1-lambdafzm)+*/U[i + DUMMY_NUM2][j - 1 + DUMMY_NUM2][k]/**lambdafzm*/;

		Urp1[k]=/*U[i][j][k]*(1-lambdafrp)+*/U[i + 1  + DUMMY_NUM2][j + DUMMY_NUM2][k]/**lambdafrp*/;

		Urm1[k]=/*U[i][j][k]*(1-lambdafrm)+*/U[i - 1  + DUMMY_NUM2][j + DUMMY_NUM2][k]/**lambdafrm*/;
	}

	for(int k = 0; k < NInitial::getNcomp(); k++)
	{
		grad[0][k]=(wyp1*pow(deltaryp1,2)*(Urp1[k]-U[i + DUMMY_NUM2][j + DUMMY_NUM2][k])-wym1*pow(deltarym1,2)*(Urm1[k]-U[i + DUMMY_NUM2][j + DUMMY_NUM2][k]))/
			(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
		grad[1][k]=(wxp1*pow(deltarxp1,2)*(Uzp1[k]-U[i + DUMMY_NUM2][j + DUMMY_NUM2][k])-wxm1*pow(deltarxm1,2)*(Uzm1[k]-U[i + DUMMY_NUM2][j + DUMMY_NUM2][k]))/
			(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));
	}

	delete []Uzp1;
	delete []Urp1;
	delete []Uzm1;
	delete []Urm1;
}
void NMethods::get_grad1Visc(const int i, const int j, const double ***U, double **grad, const char *bz, const char *br,
	const double *hr, const double *hz)
{
	double Sr1[2], Sr2[2], Sz[2], *Uzp1, *Urp1, lambdafrp, lambdafzm, lambdafrm, lambdafzp, vol,
		*Uzm1, *Urm1, *U1;

	Uzp1 = new double[NInitial::getNcomp()];
	Urp1 = new double[NInitial::getNcomp()];
	Uzm1 = new double[NInitial::getNcomp()];
	Urm1 = new double[NInitial::getNcomp()];
	U1 = new double[NInitial::getNcomp()];

	if((i>=0) && (j>=0))
	{
		Sr1[0]=hz[j];
		Sr2[0]=hz[j];
		Sz[1]=hr[i];
		vol=hz[j]*hr[i];
	}
	else if((i<0) && (j<0))
	{
		Sr1[0]=hz[abs(j)-1];
		Sr2[0]=hz[abs(j)-1];
		Sz[1]=hr[abs(i)-1];
		vol=hz[abs(j)-1]*hr[abs(i)-1];
	}
	else if((i<0) && (j>=0))
	{
		Sr1[0]=hz[j];
		Sr2[0]=hz[j];
		Sz[1]=hr[abs(i)-1];
		vol=hz[j]*hr[abs(i)-1];
	}
	else if((j<0) && (i>=0))
	{
		Sz[1]=hr[i];
		Sr1[0]=hz[abs(j)-1];
		Sr2[0]=hz[abs(j)-1];
		vol=hz[abs(j)-1]*hr[i];
	}

	Sr1[1]=0;
	Sr2[1]=0;
	Sz[0]=0;

			//double Ust[Ncomp],U1r[Ncomp],U1z[Ncomp];

				//for(int l=0;l<=Ncomp-1;l++)
				//Ust[l]=U[i][j][l];

				//if(br=="slip")
				//slip(Ust,U1r,"r");
				//else
				//periodic(i,j,const_cast<double***>(U),U1r,"r");

				//if(bz=="slip")
				//slip(Ust,U1z,"z");
				//else
				//periodic(i,j,const_cast<double***>(U),U1z,"z");



			//if((i!=0) && (i!=xmax-1) && (j!=0) && (j!=ymax-1))
			//{	
	if(i>0)
	{
		lambdafrp=hr[i]/(hr[i]+hr[i+1]);
		lambdafrm=hr[i]/(hr[i]+hr[i-1]);
	}
	else if(i==0)
	{
		lambdafrp=hr[i]/(hr[i]+hr[i+1]);
		lambdafrm=hr[i]/(hr[i]+hr[0]);
	}
	else if(i==-1)
	{
		lambdafrp=hr[abs(i)-1]/(hr[abs(i)-1]+hr[i+1]);
		lambdafrm=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i-1)-1]);
	}
	else 
	{
		lambdafrp=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i+1)-1]);
		lambdafrm=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i-1)-1]);
	}

	if(j>0)
	{
		lambdafzp=hz[j]/(hz[j]+hz[j+1]);
		lambdafzm=hz[j]/(hz[j]+hz[j-1]);
	}
	else if(j==0)
	{
		lambdafzp=hz[j]/(hz[j]+hz[j+1]);
		lambdafzm=hz[j]/(hz[j]+hz[abs(j-1)-1]);
	}
	else if(j==-1)
	{
		lambdafzp=hz[abs(j)-1]/(hz[abs(j)-1]+hz[j+1]);
		lambdafzm=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j-1)-1]);
	}
	else 
	{
		lambdafzp=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j+1)-1]);
		lambdafzm=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j-1)-1]);
	}

	double wxp1,wxm1,wyp1,wym1,deltarxp1,deltarym1,deltarxm1,deltaryp1;

	if(i>-1)
	{
		if(j>0)
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[j]+hz[j-1])*.5);
		}
		else if(j==0)
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[j]+hz[abs(j-1)-1])*.5);
		}
		else if(j==-1)
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[abs(j)-1]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
		else
		{
			wxp1=lambdafzp*hr[i]/(vol*(hz[abs(j)-1]+hz[abs(j+1)-1])*.5);
			wxm1=lambdafzm*hr[i]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
	}
	else
	{
		if(j>0)
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[j]+hz[j-1])*.5);
		}
		else if(j==0)
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[j]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[j]+hz[abs(j-1)-1])*.5);
		}
		else if(j==-1)
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[j+1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
		else
		{
			wxp1=lambdafzp*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[abs(j+1)-1])*.5);
			wxm1=lambdafzm*hr[abs(i)-1]/(vol*(hz[abs(j)-1]+hz[abs(j-1)-1])*.5);
		}
	}


	if(j>-1)
	{
		if(i>0)
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[i]+hr[i-1])*.5);
		}
		else if(i==0)
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[i]+hr[abs(i-1)-1])*.5);
		}
		else if(i==-1)
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[abs(i)-1]+hr[i+1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
		else
		{
			wyp1=lambdafrp*hz[j]/(vol*(hr[abs(i)-1]+hr[abs(i+1)-1])*.5);
			wym1=lambdafrm*hz[j]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
	}
	else
	{
		if(i>0)
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[i]+hr[i-1])*.5);
		}
		else if(i==0)
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[i]+hr[i+1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[i]+hr[abs(i-1)-1])*.5);
		}
		else if(i==-1)
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[i+1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
		else
		{
			wyp1=lambdafrp*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[abs(i+1)-1])*.5);
			wym1=lambdafrm*hz[abs(j)-1]/(vol*(hr[abs(i)-1]+hr[abs(i-1)-1])*.5);
		}
	}
				
	if(j>0)
	{
		deltarxp1=(hz[j]+hz[j+1])*.5;
		deltarxm1=(hz[j]+hz[j-1])*.5;
	}
	else if(j==0)
	{
		deltarxp1=(hz[j]+hz[j+1])*.5;
		deltarxm1=(hz[j]+hz[abs(j-1)-1])*.5;
	}
	else if(j==-1)
	{
		deltarxp1=(hz[abs(j)-1]+hz[j+1])*.5;
		deltarxm1=(hz[abs(j)-1]+hz[abs(j-1)-1])*.5;
	}
	else
	{
		deltarxp1=(hz[abs(j)-1]+hz[abs(j+1)-1])*.5;
		deltarxm1=(hz[abs(j)-1]+hz[abs(j-1)-1])*.5;
	}

	if(i>0)
	{
		deltaryp1=(hr[i]+hr[i+1])*.5;
		deltarym1=(hr[i]+hr[i-1])*.5;
	}
	else if(i==0)
	{
		deltaryp1=(hr[i]+hr[i+1])*.5;
		deltarym1=(hr[i]+hr[abs(i-1)-1])*.5;
	}
	else if(i==-1)
	{
		deltaryp1=(hr[abs(i)-1]+hr[i+1])*.5;
		deltarym1=(hr[abs(i)-1]+hr[abs(i-1)-1])*.5;
	}
	else
	{
		deltaryp1=(hr[abs(i)-1]+hr[abs(i+1)-1])*.5;
		deltarym1=(hr[abs(i)-1]+hr[abs(i-1)-1])*.5;
	}

	if((i > -1) & (j > -1) & (i < NInitial::get_xmax()) & (j < NInitial::get_ymax()))
	{
		for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
		{
			if(j != NInitial::get_ymax() - 1)
				Uzp1[k]=/*U[i][j][k]*(1-lambdafzp)+*/U[i][j+1][k]/**lambdafzp*/;
			else
			{
				double *Ust, *U1z;

				Ust = new double[NInitial::getNcomp()];
				U1z = new double[NInitial::getNcomp()];

				NBoundaryCond *bound = new NBoundaryCond;

				for(int l = 0; l <= NInitial::getNcomp() - 1; l++)
					Ust[l]=U[i][j][l];

				if(bz=="slip")
					bound->slipVisc(Ust,U1z,"z");
				else if(bz=="per")
					bound->periodic(i,j,const_cast<double***>(U),U1z,"z");
				else
					bound->inflow(i, NInitial::get_pin(), NInitial::get_roin(), NInitial::get_urin(), NInitial::get_uzin(), NInitial::get_uphiin(), U1z, hr);

				Uzp1[k]=/*U[i][j][k]*(1-lambdafzp)+*/U1z[k]/**lambdafzp*/;

				delete []Ust;
				delete []U1z;

				delete bound;
			}

			if(j!=0)
				Uzm1[k]=/*U[i][j][k]*(1-lambdafzm)+*/U[i][j-1][k]/**lambdafzm*/;
			else
			{
				double *Ust, *U1z;

				Ust = new double[NInitial::getNcomp()];
				U1z = new double[NInitial::getNcomp()];

				NBoundaryCond *bound = new NBoundaryCond;

				for(int l=0;l<=NInitial::getNcomp()-1;l++)
					Ust[l]=U[i][j][l];

				if(bz=="slip")
					bound->slipVisc(Ust,U1z,"z");
				else if(bz=="per")
					bound->periodic(i,j,const_cast<double***>(U),U1z,"z");
				else
					bound->outflow(Ust,U1z);

				Uzm1[k]=/*U[i][j][k]*(1-lambdafzm)+*/U1z[k]/**lambdafzm*/;

				delete []Ust;
				delete []U1z;

				delete bound;
			}

			if(i!=NInitial::get_xmax()-1)
				Urp1[k]=/*U[i][j][k]*(1-lambdafrp)+*/U[i+1][j][k]/**lambdafrp*/;
			else
			{
				double *Ust, *U1r;

				Ust = new double[NInitial::getNcomp()];
				U1r = new double[NInitial::getNcomp()];

				NBoundaryCond *bound = new NBoundaryCond;

				for(int l=0;l<=NInitial::getNcomp()-1;l++)
					Ust[l]=U[i][j][l];

				if(br=="slip")
					bound->slipVisc(Ust,U1r,"r");
				else
					bound->periodic(i,j,const_cast<double***>(U),U1r,"r");

				Urp1[k]=/*U[i][j][k]*(1-lambdafrp)+*/U1r[k]/**lambdafrp*/;

				delete []Ust;
				delete []U1r;

				delete bound;
			}

			if(i!=0)
				Urm1[k]=/*U[i][j][k]*(1-lambdafrm)+*/U[i-1][j][k]/**lambdafrm*/;
			else
			{
				double *Ust, *U1r;

				Ust = new double[NInitial::getNcomp()];
				U1r = new double[NInitial::getNcomp()];

				NBoundaryCond *bound = new NBoundaryCond;

				for(int l=0;l<=NInitial::getNcomp()-1;l++)
					Ust[l]=U[i][j][l];

				if(br=="slip")
					bound->slipinnerVisc(Ust,U1r);
				else
					bound->periodic(i,j,const_cast<double***>(U),U1r,"r");

				Urm1[k]=/*U[i][j][k]*(1-lambdafrm)+*/U1r[k]/**lambdafrm*/;

				delete []Ust;
				delete []U1r;

				delete bound;
			}
		}

		for(int k=0;k<NInitial::getNcomp();k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(Urp1[k]-U[i][j][k])-wym1*pow(deltarym1,2)*(Urm1[k]-U[i][j][k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(Uzp1[k]-U[i][j][k])-wxm1*pow(deltarxm1,2)*(Uzm1[k]-U[i][j][k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));
		}
	}
//center cell-----------------------------------	

	double *U1adjrp, *U1adjzp, *U1adjzm, *U1adjrm;

	U1adjrp = new double[NInitial::getNcomp()];
	U1adjzp = new double[NInitial::getNcomp()];
	U1adjzm = new double[NInitial::getNcomp()];
	U1adjrm = new double[NInitial::getNcomp()];

	NBoundaryCond *bound = new NBoundaryCond;

	if((i<0) && (j>-1) && (j<NInitial::get_ymax()))
	{
		if(br=="slip")
			bound->slipinnerVisc(U[abs(i)-1][j],U1);
		else
			bound->periodic(abs(i)-1,j,const_cast<double***>(U),U1,"r");


		if((i+1) < 0)
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i+1)-1][j],U1adjrp);
			else
				bound->periodic(abs(i+1)-1,j,const_cast<double***>(U),U1adjrp,"r");
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjrp[l]=U[i+1][j][l];


		if(br=="slip")
			bound->slipinnerVisc(U[abs(i-1)-1][j],U1adjrm);
		else
			bound->periodic(abs(i-1)-1,j,const_cast<double***>(U),U1adjrm,"r");


		if((j+1) < NInitial::get_ymax())
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][j+1],U1adjzp);
			else
				bound->periodic(abs(i)-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][2*NInitial::get_ymax()-j-2],U1adjzp);
			else
				bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");
		}


		if((j-1) > -1)
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][j-1],U1adjzm);
			else
				bound->periodic(abs(i)-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][abs(j-1)-1],U1adjzm);
			else
				bound->periodic(abs(i)-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));

						//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}
	}
	else if((i<0) && (j<0))
	{
		if(br=="slip")
			bound->slipinnerVisc(U[abs(i)-1][abs(j)-1],U1);
		else
			bound->periodic(abs(i)-1,abs(j)-1,const_cast<double***>(U),U1,"r");

		if((i+1) < 0)
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i+1)-1][abs(j)-1],U1adjrp);
			else
				bound->periodic(abs(i+1)-1,abs(j)-1,const_cast<double***>(U),U1adjrp,"r");
		}
		else
		{
			if(bz=="slip")
				bound->slipVisc(U[i+1][abs(j)-1],U1adjrp,"z");
			else if(bz=="per")
				bound->periodic(i+1,abs(j)-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->outflow(U[i+1][abs(j)-1],U1adjrp);
		}
		
		if((j+1) < 0)
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][abs(j+1)-1],U1adjzp);
			else
				bound->periodic(abs(i)-1,abs(j+1)-1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][j+1],U1adjzp);
			else
				bound->periodic(abs(i)-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}

		if(br=="slip")
			bound->slipinnerVisc(U[abs(i)-1][abs(j-1)-1],U1adjzm);
		else
			bound->periodic(abs(i)-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");

					
		if(br=="slip")
			bound->slipinnerVisc(U[abs(i-1)-1][abs(j)-1],U1adjrm);
		else
			bound->periodic(abs(i-1)-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"r");


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));

						//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}


	}
	else if((i<0) && (j>NInitial::get_ymax()-1))
	{
		if(br=="slip")
			bound->slipinnerVisc(U[abs(i)-1][2*NInitial::get_ymax()-j-1],U1);
		else
			bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1,"r");

		if((j-1) > (NInitial::get_ymax()-1))
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][2*NInitial::get_ymax()-j],U1adjzm);
			else
				bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i)-1][j-1],U1adjzm);
			else
				bound->periodic(abs(i)-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}

		if((i+1) < 0)
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i+1)-1][2*NInitial::get_ymax()-j-1],U1adjrp);
			else
				bound->periodic(abs(i+1)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"r");
		}
		else
		{
			if(bz=="slip")
				bound->slipVisc(U[i+1][2*NInitial::get_ymax()-j-1],U1adjrp,"z");
			else if(bz=="per")
				bound->periodic(i+1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrp,hr);
		}

		if(br=="slip")
			bound->slipinnerVisc(U[abs(i)-1][2*NInitial::get_ymax()-j-2],U1adjzp);
		else
			bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");
		
		if(br=="slip")
			bound->slipinnerVisc(U[abs(i-1)-1][2*NInitial::get_ymax()-j-1],U1adjrm);
		else
			bound->periodic(abs(i-1)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"r");
		
		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));

						//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}


	}
	else if((i>-1) && (i<NInitial::get_xmax()) && (j>NInitial::get_ymax()-1))
	{
		if(bz=="slip")
			bound->slipVisc(U[i][2*NInitial::get_ymax()-j-1],U1,"z");
		else if(bz=="per")
			bound->periodic(i,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1,"z");
		else
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
			NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);

		if((j-1) > (NInitial::get_ymax()-1))
		{
			if(bz=="slip")
				bound->slipVisc(U[i][2*NInitial::get_ymax()-j],U1adjzm,"z");
			else if(bz=="per")
				bound->periodic(i,2*NInitial::get_ymax()-j,const_cast<double***>(U),U1adjzm,"z");
			else
				bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1adjzm,hr);
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjzm[l]=U[i][j-1][l];

		if((i+1) < NInitial::get_xmax())
		{
			if(bz=="slip")
				bound->slipVisc(U[i+1][2*NInitial::get_ymax()-j-1],U1adjrp,"z");
			else if(bz=="per")
				bound->periodic(i+1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrp,hr);
		}
		else
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-2][2*NInitial::get_ymax()-j-1],U1adjrp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-2,2*NInitial::get_ymax()-j-1,
				const_cast<double***>(U),U1adjrp,"r");
		}


		if((i-1) > -1)
		{
			if(bz=="slip")
				bound->slipVisc(U[i-1][2*NInitial::get_ymax()-j-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrm,hr);
		}
		else
		{
			if(br=="slip")
				bound->slipVisc(U[abs(i-1)-1][2*NInitial::get_ymax()-j-1],U1adjrm,"r");
			else
				bound->periodic(abs(i-1)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"r");
		}


		if(bz=="slip")
			bound->slipVisc(U[i][2*NInitial::get_ymax()-j-2],U1adjzp,"z");
		else if(bz=="per")
			bound->periodic(i,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"z");
		else
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
			NInitial::get_uzin(),NInitial::get_uphiin(),U1adjzp,hr);



		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));
			
			//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}


	}
	else if((i>NInitial::get_xmax()-1) && (j>NInitial::get_ymax()-1))
	{
		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j-1],U1,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j-1,
			const_cast<double***>(U),U1,"r");
					
		if((j-1) > (NInitial::get_ymax()-1))
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][j-1],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}

		if((i-1) > (NInitial::get_xmax()-1))
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-2][2*NInitial::get_ymax()-j-1],U1adjrm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-2,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"r");
		}
		else
		{
			if(bz=="slip")
				bound->slipVisc(U[i-1][2*NInitial::get_ymax()-j-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrm,hr);
		}

		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j-2],U1adjzp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");

		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-2][2*NInitial::get_ymax()-j-1],U1adjrp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-2,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));

						//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}


	}
	else if((i>NInitial::get_xmax()-1) && (j<NInitial::get_ymax()) && (j>-1))
	{
		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-1][j],U1,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,j,const_cast<double***>(U),U1,"r");

		if((j+1) < NInitial::get_ymax())
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][j+1],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j-2],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");
		}

		if((j-1) > -1)
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][j-1],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][abs(j-1)-1],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");
		}

		if((i-1) > (NInitial::get_xmax()-1))
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i][j],U1adjrm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i,j,const_cast<double***>(U),U1adjrm,"r");
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjrm[l]=U[NInitial::get_xmax()-1][j][l];


		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-2][j],U1adjrp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-2,j,const_cast<double***>(U),U1adjrp,"r");



		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));
			
			//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}


	}
	else if((i>NInitial::get_xmax()-1) && (j<0))
	{
		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-1][abs(j)-1],U1,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,abs(j)-1,const_cast<double***>(U),U1,"r");

		if((j+1) < 0)
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][abs(j+1)-1],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,abs(j+1)-1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-1][j+1],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}

					
		if((i-1) > (NInitial::get_xmax()-1))
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i][abs(j)-1],U1adjrm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i,abs(j)-1,const_cast<double***>(U),U1adjrm,"r");
		}
		else
		{
			if(bz=="slip")
				bound->slipVisc(U[i-1][abs(j)-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->outflow(U[0][abs(j)-1],U1adjrm);
		}


		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-2][abs(j)-1],U1adjrp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-2,abs(j)-1,const_cast<double***>(U),U1adjrp,"r");


		if(br=="slip")
			bound->slipVisc(U[2*NInitial::get_xmax()-i-1][abs(j-1)-1],U1adjzm,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));

			//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}


	}
	else if((i<NInitial::get_xmax()) && (i>-1) && (j<0))
	{
		if(bz=="slip")
			bound->slipVisc(U[i][abs(j)-1],U1,"z");
		else if(bz=="per")
			bound->periodic(i,abs(j)-1,const_cast<double***>(U),U1,"z");
		else
			bound->outflow(U[i][0],U1);


		if((j+1) < 0)
		{
			if(bz=="slip")
				bound->slipVisc(U[i][abs(j+1)-1],U1adjzp,"z");
			else if(bz=="per")
				bound->periodic(i,abs(j+1)-1,const_cast<double***>(U),U1adjzp,"z");
			else
				bound->outflow(U[i][0],U1adjzp);
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjzp[l]=U[i][j+1][l];


		if((i-1) > -1)
		{
			if(bz=="slip")
				bound->slipVisc(U[i-1][abs(j)-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->outflow(U[i-1][0],U1adjrm);
		}
		else
		{
			if(br=="slip")
				bound->slipinnerVisc(U[abs(i-1)-1][abs(j)-1],U1adjrm);
			else
				bound->periodic(abs(i-1)-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"r");
		}


		if((i+1) < NInitial::get_xmax())
		{
			if(bz=="slip")
				bound->slipVisc(U[i+1][abs(j)-1],U1adjrp,"z");
			else if(bz=="per")
				bound->periodic(i+1,abs(j)-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->outflow(U[i+1][0],U1adjrp);
		}
		else
		{
			if(br=="slip")
				bound->slipVisc(U[2*NInitial::get_xmax()-i-2][abs(j)-1],U1adjrp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-2,abs(j)-1,const_cast<double***>(U),U1adjrp,"r");
		}	


		if(bz=="slip")
			bound->slipVisc(U[i][abs(j-1)-1],U1adjzm,"z");
		else if(bz=="per")
			bound->periodic(i,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"z");
		else
			bound->outflow(U[i][0],U1adjzm);


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			grad[0][k]=(wyp1*pow(deltaryp1,2)*(U1adjrp[k]-U1[k])-wym1*pow(deltarym1,2)*(U1adjrm[k]-U1[k]))/
				(wyp1*pow(deltaryp1,2)+wym1*pow(deltarym1,2));
			grad[1][k]=(wxp1*pow(deltarxp1,2)*(U1adjzp[k]-U1[k])-wxm1*pow(deltarxm1,2)*(U1adjzm[k]-U1[k]))/
				(wxp1*pow(deltarxp1,2)+wxm1*pow(deltarxm1,2));

						//Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp;
						//Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm;
						//Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp;
						//Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm;
		}

		//for(int l=0;l<Ncomp;l++)
				//{
				//	if(Uzp1[l]!=Uzp1[l])
				//		cout<<"zp"<<' ';
				//	if(Uzm1[l]!=Uzm1[l])
				//		cout<<"zm"<<' ';
				//	if(Urp1[l]!=Urp1[l])
				//		cout<<"rp"<<' ';
				//	if(Urm1[l]!=Urm1[l])
				//		cout<<"rm"<<' ';

				//}


	}

	delete []Uzp1;
	delete []Urp1;
	delete []Uzm1;
	delete []Urm1;
	delete []U1;

	delete []U1adjrp;
	delete []U1adjzp;
	delete []U1adjzm;
	delete []U1adjrm;

	delete bound;
}
void NMethods::get_grad2(const int i, const int j, const double ***U, double **grad, const char *bz, const char *br,
	const double *hr, const double *hz)
{		
	double Sr1[2],Sr2[2],Sz[2],*Uzp1,*Urp1,*Urm1,*Uzm1,lambdafrp,lambdafrm,lambdafzp,lambdafzm,vol;
	
	Uzp1 = new double[NInitial::getNcomp()];
	Urp1 = new double[NInitial::getNcomp()];
	Urm1 = new double[NInitial::getNcomp()];
	Uzm1 = new double[NInitial::getNcomp()];

	if(i>0)
	{
		lambdafrp=hr[i]/(hr[i]+hr[i+1]);
		lambdafrm=hr[i]/(hr[i]+hr[i-1]);
	}
	else if(i==0)
	{
		lambdafrp=hr[i]/(hr[i]+hr[i+1]);
		lambdafrm=hr[i]/(hr[i]+hr[0]);
	}
	else if(i==-1)
	{
		lambdafrp=hr[abs(i)-1]/(hr[abs(i)-1]+hr[i+1]);
		lambdafrm=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i-1)-1]);
	}
	else 
	{
		lambdafrp=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i+1)-1]);
		lambdafrm=hr[abs(i)-1]/(hr[abs(i)-1]+hr[abs(i-1)-1]);
	}

	if(j>0)
	{
		lambdafzp=hz[j]/(hz[j]+hz[j+1]);
		lambdafzm=hz[j]/(hz[j]+hz[j-1]);
	}
	else if(j==0)
	{
		lambdafzp=hz[j]/(hz[j]+hz[j+1]);
		lambdafzm=hz[j]/(hz[j]+hz[abs(j-1)-1]);
	}
	else if(j==-1)
	{
		lambdafzp=hz[abs(j)-1]/(hz[abs(j)-1]+hz[j+1]);
		lambdafzm=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j-1)-1]);
	}
	else 
	{
		lambdafzp=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j+1)-1]);
		lambdafzm=hz[abs(j)-1]/(hz[abs(j)-1]+hz[abs(j-1)-1]);
	}

	if((i>=0) && (j>=0))
	{
		Sr1[0]=hz[j];
		Sr2[0]=hz[j];
		Sz[1]=hr[i];
		vol=hz[j]*hr[i];
	}
	else if((i<0) && (j<0))
	{
		Sr1[0]=hz[abs(j)-1];
		Sr2[0]=hz[abs(j)-1];
		Sz[1]=hr[abs(i)-1];
		vol=hz[abs(j)-1]*hr[abs(i)-1];
	}
	else if((i<0) && (j>=0))
	{
		Sr1[0]=hz[j];
		Sr2[0]=hz[j];
		Sz[1]=hr[abs(i)-1];
		vol=hz[j]*hr[abs(i)-1];
	}
	else if((j<0) && (i>=0))
	{
		Sz[1]=hr[i];
		Sr1[0]=hz[abs(j)-1];
		Sr2[0]=hz[abs(j)-1];
		vol=hz[abs(j)-1]*hr[i];
	}
	
	Sr1[1]=0;
	Sr2[1]=0;
	Sz[0]=0;

	double **grU, **grUadj, **grUzp1, **grUzm1, **grUrp1, **grUrm1;

	grU = new double *[2];
	for(int i = 0; i < 2; i++)
		grU[i] = new double[NInitial::getNcomp()];
	grUadj = new double* [2];
	for(int i = 0; i < 2; i++)
		grUadj[i] = new double[NInitial::getNcomp()];
	grUzp1 = new double* [2];
	for(int i = 0; i < 2; i++)
		grUzp1[i] = new double[NInitial::getNcomp()];
	grUzm1 = new double* [2];
	for(int i = 0; i < 2; i++)
		grUzm1[i] = new double[NInitial::getNcomp()];
	grUrp1 = new double* [2];
	for(int i = 0; i < 2; i++)
		grUrp1[i] = new double[NInitial::getNcomp()];
	grUrm1 = new double* [2];
	for(int i = 0; i < 2; i++)
		grUrm1[i] = new double[NInitial::getNcomp()];



	get_grad1(i,j,const_cast<double***>(U),grU,bz,br,hr,hz);
	get_grad1(i,j+1,const_cast<double***>(U),grUadj,bz,br,hr,hz);
	
	for(int l=0;l<=NInitial::getNcomp()-1;l++)
		for(int k=0;k<=1;k++)
			grUzp1[k][l]=grU[k][l]*(1-lambdafzp)+grUadj[k][l]*lambdafzp;
					
			//for(int k=0;k<=1;k++)
			//	for(int l=0;l<Ncomp;l++)
			//	{
			//		if(grU[k][l]!=grU[k][l])
			//			cout<<"U"<<' '<<j<<' '<<i<<' ';
			//		if(grUadj[k][l]!=grUadj[k][l])
			//			cout<<"adj"<<' ';
			//	}

	get_grad1(i,j,const_cast<double***>(U),grU,bz,br,hr,hz);
	get_grad1(i,j-1,const_cast<double***>(U),grUadj,bz,br,hr,hz);
	
	for(int l=0;l<=NInitial::getNcomp()-1;l++)
		for(int k=0;k<=1;k++)
			grUzm1[k][l]=grU[k][l]*(1-lambdafzm)+grUadj[k][l]*lambdafzm;

	get_grad1(i,j,const_cast<double***>(U),grU,bz,br,hr,hz);
	get_grad1(i+1,j,const_cast<double***>(U),grUadj,bz,br,hr,hz);
	
	for(int l=0;l<=NInitial::getNcomp()-1;l++)
		for(int k=0;k<=1;k++)
			grUrp1[k][l]=grU[k][l]*(1-lambdafrp)+grUadj[k][l]*lambdafrp;

	get_grad1(i,j,const_cast<double***>(U),grU,bz,br,hr,hz);
	get_grad1(i-1,j,const_cast<double***>(U),grUadj,bz,br,hr,hz);
	
	for(int l=0;l<=NInitial::getNcomp()-1;l++)
		for(int k=0;k<=1;k++)
			grUrm1[k][l]=grU[k][l]*(1-lambdafrm)+grUadj[k][l]*lambdafrm;


	if((i>-1) && (j>-1) && (i<NInitial::get_xmax()) && (j<NInitial::get_ymax()))
	{
		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			if(j!=NInitial::get_ymax()-1)
				Uzp1[k]=U[i][j][k]*(1-lambdafzp)+U[i][j+1][k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
			else
			{
				double *Ust, *U1z;

				Ust = new double[NInitial::getNcomp()];
				U1z = new double[NInitial::getNcomp()];

				NMgdBoundaryCond *bound = new NMgdBoundaryCond;

				for(int l=0;l<=NInitial::getNcomp()-1;l++)
					Ust[l]=U[i][j][l];

				if(bz=="slip")
					bound->slip(Ust,U1z,"z");
				else if(bz=="per")
					bound->periodic(i,j,const_cast<double***>(U),U1z,"z");
				else
					bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
					NInitial::get_uzin(),NInitial::get_uphiin(),U1z,hr);


				Uzp1[k]=U[i][j][k]*(1-lambdafzp)+U1z[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);

				delete []Ust;
				delete []U1z;
				delete bound;
			}

			if(j!=0)
				Uzm1[k]=U[i][j][k]*(1-lambdafzm)+U[i][j-1][k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			else
			{
				double *Ust, *U1z;

				Ust = new double[NInitial::getNcomp()];
				U1z = new double[NInitial::getNcomp()];

				NMgdBoundaryCond *bound = new NMgdBoundaryCond;

				for(int l=0;l<=NInitial::getNcomp()-1;l++)
					Ust[l]=U[i][j][l];

				if(bz=="slip")
					bound->slip(Ust,U1z,"z");
				else if(bz=="per")
					bound->periodic(i,j,const_cast<double***>(U),U1z,"z");
				else
					bound->outflow(Ust,U1z);


				Uzm1[k]=U[i][j][k]*(1-lambdafzm)+U1z[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[0]*.5*lambdafzm-hz[j]*.5);

				delete []Ust;
				delete []U1z;
				delete bound;
			}

			if(i!=NInitial::get_xmax()-1)
				Urp1[k]=U[i][j][k]*(1-lambdafrp)+U[i+1][j][k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
			else
			{
				double *Ust, *U1r;

				Ust = new double[NInitial::getNcomp()];
				U1r = new double[NInitial::getNcomp()];

				NMgdBoundaryCond *bound = new NMgdBoundaryCond;

				for(int l=0;l<=NInitial::getNcomp()-1;l++)
					Ust[l]=U[i][j][l];

				if(br=="slip")
					bound->slip(Ust,U1r,"r");
				else
					bound->periodic(i,j,const_cast<double***>(U),U1r,"r");


				Urp1[k]=U[i][j][k]*(1-lambdafrp)+U1r[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);

				delete []Ust;
				delete []U1r;
				delete bound;
			}

			if(i!=0)
				Urm1[k]=U[i][j][k]*(1-lambdafrm)+U[i-1][j][k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			else
			{
				double *Ust, *U1r;

				Ust = new double[NInitial::getNcomp()];
				U1r = new double[NInitial::getNcomp()];

				NMgdBoundaryCond *bound = new NMgdBoundaryCond;

				for(int l=0;l<=NInitial::getNcomp()-1;l++)
					Ust[l]=U[i][j][l];

				if(br=="slip")
					bound->slip(Ust,U1r,"r");
				else
					bound->periodic(i,j,const_cast<double***>(U),U1r,"r");


				Urm1[k]=U[i][j][k]*(1-lambdafrm)+U1r[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[0]*.5*lambdafrm-hr[i]*.5);

				delete []Ust;
				delete []U1r;
				delete bound;
			}
		}
	}

				
	double *U1, *U1adjrp, *U1adjrm, *U1adjzp, *U1adjzm;	
	
	U1 = new double[NInitial::getNcomp()];
	U1adjrp = new double[NInitial::getNcomp()];
	U1adjrm = new double[NInitial::getNcomp()];
	U1adjzp = new double[NInitial::getNcomp()];
	U1adjzm = new double[NInitial::getNcomp()];

	NMgdBoundaryCond *bound = new NMgdBoundaryCond;

	if((i<0) && (j>-1) && (j<NInitial::get_ymax()))
	{
		if(br=="slip")
			bound->slip(U[abs(i)-1][j],U1,"r");
		else
			bound->periodic(abs(i)-1,j,const_cast<double***>(U),U1,"r");

		if((i+1) < 0)
		{
			if(br=="slip")
				bound->slip(U[abs(i+1)-1][j],U1adjrp,"r");
			else
				bound->periodic(abs(i+1)-1,j,const_cast<double***>(U),U1adjrp,"r");
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjrp[l]=U[i+1][j][l];


		if(br=="slip")
			bound->slip(U[abs(i-1)-1][j],U1adjrm,"r");
		else
			bound->periodic(abs(i-1)-1,j,const_cast<double***>(U),U1adjrm,"r");


		if((j+1) < NInitial::get_ymax())
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][j+1],U1adjzp,"r");
			else
				bound->periodic(abs(i)-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][2*NInitial::get_ymax()-j-2],U1adjzp,"r");
			else
				bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");
		}


		if((j-1) > -1)
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][j-1],U1adjzm,"r");
			else
				bound->periodic(abs(i)-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][abs(j-1)-1],U1adjzm,"r");
			else
				bound->periodic(abs(i)-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");
		}
		

		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			

			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}
	}
	else if((i<0) && (j<0))
	{
		if(br=="slip")
			bound->slip(U[abs(i)-1][abs(j)-1],U1,"r");
		else
			bound->periodic(abs(i)-1,abs(j)-1,const_cast<double***>(U),U1,"r");

		if((i+1) < 0)
		{
			if(br=="slip")
				bound->slip(U[abs(i+1)-1][abs(j)-1],U1adjrp,"r");
			else
				bound->periodic(abs(i+1)-1,abs(j)-1,const_cast<double***>(U),U1adjrp,"r");
		}
		else
		{
			if(bz=="slip")
				bound->slip(U[i+1][abs(j)-1],U1adjrp,"z");
			else if(bz=="per")
				bound->periodic(i+1,abs(j)-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->outflow(U[i+1][abs(j)-1],U1adjrp);
		}


		if((j+1) < 0)
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][abs(j+1)-1],U1adjzp,"r");
			else
				bound->periodic(abs(i)-1,abs(j+1)-1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][j+1],U1adjzp,"r");
			else
				bound->periodic(abs(i)-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}


		if(br=="slip")
			bound->slip(U[abs(i)-1][abs(j-1)-1],U1adjzm,"r");
		else
			bound->periodic(abs(i)-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");
							
		if(br=="slip")
			bound->slip(U[abs(i-1)-1][abs(j)-1],U1adjrm,"r");
		else
			bound->periodic(abs(i-1)-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"r");



		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			

			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}

	}
	else if((i<0) && (j>NInitial::get_ymax()-1))
	{
		if(br=="slip")
			bound->slip(U[abs(i)-1][2*NInitial::get_ymax()-j-1],U1,"r");
		else
			bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1,"r");

		if((j-1) > (NInitial::get_ymax()-1))
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][2*NInitial::get_ymax()-j],U1adjzm,"r");
			else
				bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[abs(i)-1][j-1],U1adjzm,"r");
			else
				bound->periodic(abs(i)-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}
		
		if((i+1) < 0)
		{
			if(br=="slip")
				bound->slip(U[abs(i+1)-1][2*NInitial::get_ymax()-j-1],U1adjrp,"r");
			else
				bound->periodic(abs(i+1)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[i+1][2*NInitial::get_ymax()-j-1],U1adjrp,"z");
			else if(br=="per")
				bound->periodic(i+1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
					NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrp,hr);
		}

		if(br=="slip")
			bound->slip(U[abs(i)-1][2*NInitial::get_ymax()-j-2],U1adjzp,"r");
		else
			bound->periodic(abs(i)-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");

		if(br=="slip")
			bound->slip(U[abs(i-1)-1][2*NInitial::get_ymax()-j-1],U1adjrm,"r");
		else
			bound->periodic(abs(i-1)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"r");


		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}

			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}

	}
	else if((i>-1) && (i<NInitial::get_xmax()) && (j>NInitial::get_ymax()-1))
	{
		if(bz=="slip")
			bound->slip(U[i][2*NInitial::get_ymax()-j-1],U1,"z");
		else if(bz=="per")
			bound->periodic(i,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1,"z");
		else
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
		
		if((j-1) > (NInitial::get_ymax()-1))
		{
			if(bz=="slip")
				bound->slip(U[i][2*NInitial::get_ymax()-j],U1adjzm,"z");
			else if(bz=="per")
				bound->periodic(i,2*NInitial::get_ymax()-j,const_cast<double***>(U),U1adjzm,"z");
			else
				bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
					NInitial::get_uzin(),NInitial::get_uphiin(),U1adjzm,hr);
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjzm[l]=U[i][j-1][l];
		
		if((i+1) < NInitial::get_xmax())
		{
			if(bz=="slip")
				bound->slip(U[i+1][2*NInitial::get_ymax()-j-1],U1adjrp,"z");
			else if(bz=="per")
				bound->periodic(i+1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
					NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrp,hr);
		}
		else
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-2][2*NInitial::get_ymax()-j-1],U1adjrp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-2,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"r");
		}
		
		if((i-1) > -1)
		{
			if(bz=="slip")
				bound->slip(U[i-1][2*NInitial::get_ymax()-j-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
					NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrm,hr);
		}
		else
		{
			if(br=="slip")
				bound->slip(U[abs(i-1)-1][2*NInitial::get_ymax()-j-1],U1adjrm,"r");
			else
				bound->periodic(abs(i-1)-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"r");
		}
		
		if(bz=="slip")
			bound->slip(U[i][2*NInitial::get_ymax()-j-2],U1adjzp,"z");
		else if(bz=="per")
			bound->periodic(i,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"z");
		else
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1adjzp,hr);
		
		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			
			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}

	}
	else if((i>NInitial::get_xmax()-1) && (j>NInitial::get_ymax()-1))
	{
		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j-1],U1,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1,"r");
					
		if((j-1) > (NInitial::get_ymax()-1))
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][j-1],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}
		
		if((i-1) > (NInitial::get_xmax()-1))
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i][2*NInitial::get_ymax()-j-1],U1adjrm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"r");
		}
		else
		{
			if(bz=="slip")
				bound->slip(U[i-1][2*NInitial::get_ymax()-j-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1adjrm,hr);
		}

		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j-2],U1adjzp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");

		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-2][2*NInitial::get_ymax()-j-1],U1adjrp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-2,2*NInitial::get_ymax()-j-1,const_cast<double***>(U),U1adjrp,"r");

		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			
			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}

	}
	else if((i>NInitial::get_xmax()-1) && (j<NInitial::get_ymax()) && (j>-1))
	{
		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-1][j],U1,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,j,const_cast<double***>(U),U1,"r");

		if((j+1) < NInitial::get_ymax())
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][j+1],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][2*NInitial::get_ymax()-j-2],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,2*NInitial::get_ymax()-j-2,const_cast<double***>(U),U1adjzp,"r");
		}
		
		if((j-1) > -1)
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][j-1],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j-1,const_cast<double***>(U),U1adjzm,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][abs(j-1)-1],U1adjzm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");
		}

		if((i-1) > (NInitial::get_xmax()-1))
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i][j],U1adjrm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i,j,const_cast<double***>(U),U1adjrm,"r");
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjrm[l]=U[NInitial::get_xmax()-1][j][l];
		
		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-2][j],U1adjrp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-2,j,const_cast<double***>(U),U1adjrp,"r");

		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}


			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}

	}
	else if((i>NInitial::get_xmax()-1) && (j<0))
	{
		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-1][abs(j)-1],U1,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,abs(j)-1,const_cast<double***>(U),U1,"r");

		if((j+1) < 0)
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][abs(j+1)-1],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,abs(j+1)-1,const_cast<double***>(U),U1adjzp,"r");
		}
		else
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-1][j+1],U1adjzp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-1,j+1,const_cast<double***>(U),U1adjzp,"r");
		}

		if((i-1) > (NInitial::get_xmax()-1))
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i][abs(j)-1],U1adjrm,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i,abs(j)-1,const_cast<double***>(U),U1adjrm,"r");
		}
		else
		{
			if(bz=="slip")
				bound->slip(U[i-1][abs(j)-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->outflow(U[0][abs(j)-1],U1adjrm);
		}

		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-2][abs(j)-1],U1adjrp,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-2,abs(j)-1,const_cast<double***>(U),U1adjrp,"r");

		if(br=="slip")
			bound->slip(U[2*NInitial::get_xmax()-i-1][abs(j-1)-1],U1adjzm,"r");
		else
			bound->periodic(2*NInitial::get_xmax()-i-1,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"r");


		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}

			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}

	}
	else if((i<NInitial::get_xmax()) && (i>-1) && (j<0))
	{
		if(bz=="slip")
			bound->slip(U[i][abs(j)-1],U1,"z");
		else if(bz=="per")
			bound->periodic(i,abs(j)-1,const_cast<double***>(U),U1,"z");
		else
			bound->outflow(U[i][0],U1);


		if((j+1) < 0)
		{
			if(bz=="slip")
				bound->slip(U[i][abs(j+1)-1],U1adjzp,"z");
			else if(bz=="per")
				bound->periodic(i,abs(j+1)-1,const_cast<double***>(U),U1adjzp,"z");
			else
				bound->outflow(U[i][0],U1adjzp);
		}
		else
			for(int l=0;l<NInitial::getNcomp();l++)
				U1adjzp[l]=U[i][0][l];
		
		if((i-1) > -1)
		{
			if(bz=="slip")
				bound->slip(U[i-1][abs(j)-1],U1adjrm,"z");
			else if(bz=="per")
				bound->periodic(i-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"z");
			else
				bound->outflow(U[i-1][0],U1adjrm);
		}
		else
		{
			if(br=="slip")
				bound->slip(U[abs(i-1)-1][abs(j)-1],U1adjrm,"r");
			else
				bound->periodic(abs(i-1)-1,abs(j)-1,const_cast<double***>(U),U1adjrm,"r");
		}

		if((i+1) < NInitial::get_xmax())
		{
			if(bz=="slip")
				bound->slip(U[i+1][abs(j)-1],U1adjrp,"z");
			else if(bz=="per")
				bound->periodic(i+1,abs(j)-1,const_cast<double***>(U),U1adjrp,"z");
			else
				bound->outflow(U[i+1][0],U1adjrp);
		}
		else
		{
			if(br=="slip")
				bound->slip(U[2*NInitial::get_xmax()-i-2][abs(j)-1],U1adjrp,"r");
			else
				bound->periodic(2*NInitial::get_xmax()-i-2,abs(j)-1,const_cast<double***>(U),U1adjrp,"r");
		}	
		
		if(bz=="slip")
			bound->slip(U[i][abs(j-1)-1],U1adjzm,"z");
		else if(bz=="per")
			bound->periodic(i,abs(j-1)-1,const_cast<double***>(U),U1adjzm,"z");
		else
			bound->outflow(U[i][0],U1adjzm);

		for(int l=0;l<NInitial::getNcomp();l++)
		{
			if(U1adjzp[l]!=U1adjzp[l])
				cout<<"zp2"<<' ';
			if(U1adjzm[l]!=U1adjzm[l])
				cout<<"zm2"<<' ';
			if(U1adjrp[l]!=U1adjrp[l])
				cout<<"rp2"<<' ';
			if(U1adjrm[l]!=U1adjrm[l])
				cout<<"rm2"<<' ';
			if(U1[l]!=U1[l])
				cout<<"12"<<' ';
		}

		for(int k=0;k<NInitial::getNcomp();k++)
		{
			if(j>0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[j-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==0)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[j]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[j]+hz[abs(j-1)-1]*.5*lambdafzm-hz[j]*.5);
			}
			else if(j==-1)
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[j+1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}
			else 
			{
				Uzp1[k]=U1[k]*(1-lambdafzp)+U1adjzp[k]*lambdafzp+grUzp1[1][k]*(
					.5*hz[abs(j)-1]*(1-lambdafzp)-hz[abs(j+1)-1]*.5*lambdafzp);
				Uzm1[k]=U1[k]*(1-lambdafzm)+U1adjzm[k]*lambdafzm-grUzm1[1][k]*(
					.5*lambdafzm*hz[abs(j)-1]+hz[abs(j-1)-1]*.5*lambdafzm-hz[abs(j)-1]*.5);
			}

			if(i>0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[i-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==0)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[i]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[i]+hr[abs(i-1)-1]*.5*lambdafrm-hr[i]*.5);
			}
			else if(i==-1)
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[i+1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
			else 
			{
				Urp1[k]=U1[k]*(1-lambdafrp)+U1adjrp[k]*lambdafrp+grUrp1[0][k]*(
					.5*hr[abs(i)-1]*(1-lambdafrp)-hr[abs(i+1)-1]*.5*lambdafrp);
				Urm1[k]=U1[k]*(1-lambdafrm)+U1adjrm[k]*lambdafrm-grUrm1[0][k]*(
					.5*lambdafrm*hr[abs(i)-1]+hr[abs(i-1)-1]*.5*lambdafrm-hr[abs(i)-1]*.5);
			}
		}

	}



	for(int k=0;k<=NInitial::getNcomp()-1;k++)
		for(int l=0;l<=1;l++)
			grad[l][k]=1/vol*(Urm1[k]*Sr1[l]+Urp1[k]*Sr2[l]+Uzm1[k]*Sz[l]+Uzp1[k]*Sz[l]);

	delete bound;
	delete []U1;
	delete []U1adjrp;
	delete []U1adjrm;
	delete []U1adjzp;
	delete []U1adjzm;

	for(int i = 0; i < 2; i++)
		delete []grU[i];
	delete []grU;

	for(int i = 0; i < 2; i++)
		delete []grUadj[i];
	delete []grUadj;

	for(int i = 0; i < 2; i++)
		delete []grUzp1[i];
	delete []grUzp1;

	for(int i = 0; i < 2; i++)
		delete []grUzm1[i];
	delete []grUzm1;

	for(int i = 0; i < 2; i++)
		delete []grUrp1[i];
	delete []grUrp1;

	for(int i = 0; i < 2; i++)
		delete []grUrm1[i];
	delete []grUrm1;


	delete []Uzp1;
	delete []Urp1;
	delete []Urm1;
	delete []Uzm1;

}
void NMethods::getGradH(const int i, const int j, const double ***U, const double *hr,
	const double *hz, const char *bz, const char *br, double gradH[2])
{
	double **gradUij, gradVelosity[2][3], E, grE[2], gre[2], *U1;

	gradUij = new double *[2];
	for(int i = 0; i < 2; i++)
		gradUij[i] = new double[NInitial::getNcomp()];
	U1 = new double [NInitial::getNcomp()];

	NMgdBoundaryCond *bound = new NMgdBoundaryCond;

	if(i>-1 && j>-1 && i<NInitial::get_xmax() && j<NInitial::get_ymax())
	{
		get_grad1Visc(i,j,U,gradUij,bz,br,hr,hz);
		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U[i][j][0]*(gradUij[k][l+1]-gradUij[k][0]*U[i][j][l+1]/U[i][j][0]);//k - компонента градиента (r,z)
				//l - компонента скорости (r,z,phi)
		E=press(U[i][j])/(NInitial::get_sigma()*U[i][j][0])+.5*(pow(U[i][j][1]/U[i][j][0],2)+
			pow(U[i][j][2]/U[i][j][0],2)+pow(U[i][j][3]/U[i][j][0],2));

		for(int k=0;k<2;k++)
		{
			grE[k]=(gradUij[k][4]-E*gradUij[k][0])/U[i][j][0];
			gre[k]=grE[k]-U[i][j][1]/U[i][j][0]*gradVelosity[k][0]-U[i][j][2]/U[i][j][0]*gradVelosity[k][1]-
				U[i][j][3]/U[i][j][0]*gradVelosity[k][2];
			gradH[k]=(NInitial::get_sigma() + 1)*gre[k];
		}
	}
	else if(i==-1)
	{
		if(br=="slip")
			bound->slipinnerVisc(U[0][j],U1);
		else
			{/*r-conditions*/}

		get_grad1Visc(i,j,U,gradUij,bz,br,hr,hz);
		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента градиента (r,z)
				//l - компонента скорости (r,z,phi)
		E=press(U1)/(NInitial::get_sigma()*U1[0])+.5*(pow(U1[1]/U1[0],2)+
			pow(U1[2]/U1[0],2)+pow(U1[3]/U1[0],2));

		for(int k=0;k<2;k++)
		{
			grE[k]=(gradUij[k][4]-E*gradUij[k][0])/U1[0];
			gre[k]=grE[k]-U1[1]/U1[0]*gradVelosity[k][0]-U1[2]/U1[0]*gradVelosity[k][1]-
			U1[3]/U1[0]*gradVelosity[k][2];
			gradH[k]=(NInitial::get_sigma() + 1)*gre[k];
		}
	}
	else if(i==NInitial::get_xmax())
	{
		if(br=="slip")
			bound->slipVisc(U[NInitial::get_xmax()-1][j],U1,"r");
		else
			{/*r-conditions*/}

		get_grad1Visc(i,j,U,gradUij,bz,br,hr,hz);
		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента градиента (r,z)
				//l - компонента скорости (r,z,phi)
		E=press(U1)/(NInitial::get_sigma()*U1[0])+.5*(pow(U1[1]/U1[0],2)+
			pow(U1[2]/U1[0],2)+pow(U1[3]/U1[0],2));

		for(int k=0;k<2;k++)
		{
			grE[k]=(gradUij[k][4]-E*gradUij[k][0])/U1[0];
			gre[k]=grE[k]-U1[1]/U1[0]*gradVelosity[k][0]-U1[2]/U1[0]*gradVelosity[k][1]-
			U1[3]/U1[0]*gradVelosity[k][2];
			gradH[k]=(NInitial::get_sigma() + 1)*gre[k];
		}
	}
	else if(j==-1)
	{
		if(bz=="slip")
			bound->slipVisc(U[i][0],U1,"z");
		else if(bz=="per")
			bound->periodic(i,0,const_cast<double***>(U),U1,"z");
		else 
			bound->outflow(U[i][0],U1);

		get_grad1Visc(i,j,U,gradUij,bz,br,hr,hz);
		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента градиента (r,z)
				//l - компонента скорости (r,z,phi)
		E=press(U1)/(NInitial::get_sigma()*U1[0])+.5*(pow(U1[1]/U1[0],2)+
			pow(U1[2]/U1[0],2)+pow(U1[3]/U1[0],2));

		for(int k=0;k<2;k++)
		{
			grE[k]=(gradUij[k][4]-E*gradUij[k][0])/U1[0];
			gre[k]=grE[k]-U1[1]/U1[0]*gradVelosity[k][0]-U1[2]/U1[0]*gradVelosity[k][1]-
			U1[3]/U1[0]*gradVelosity[k][2];
			gradH[k]=(NInitial::get_sigma() + 1)*gre[k];
		}
	}
	else if(j==NInitial::get_ymax())
	{
		if(bz=="slip")
			bound->slipVisc(U[i][NInitial::get_ymax()-1],U1,"z");
		else if(bz=="per")
			bound->periodic(i,NInitial::get_ymax()-1,const_cast<double***>(U),U1,"z");
		else 
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
			NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);

		get_grad1Visc(i,j,U,gradUij,bz,br,hr,hz);
		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента градиента (r,z)
				//l - компонента скорости (r,z,phi)
		E=press(U1)/(NInitial::get_sigma()*U1[0])+.5*(pow(U1[1]/U1[0],2)+
			pow(U1[2]/U1[0],2)+pow(U1[3]/U1[0],2));

		for(int k=0;k<2;k++)
		{
			grE[k]=(gradUij[k][4]-E*gradUij[k][0])/U1[0];
			gre[k]=grE[k]-U1[1]/U1[0]*gradVelosity[k][0]-U1[2]/U1[0]*gradVelosity[k][1]-
			U1[3]/U1[0]*gradVelosity[k][2];
			gradH[k]=(NInitial::get_sigma() + 1)*gre[k];
		}
	}

	for(int i = 0; i < 2; i++)
		delete []gradUij[i];
	delete []gradUij;
	delete []U1;

	delete bound;
}
void NMethods::SourseVisc(const double *U, double *S, double (&tau)[3][3], const double Hr)
{
	double Pr;

	Pr=NInitial::getCp()*NInitial::get_mu()/NInitial::get_kappa();
	S[0]=0;
	S[1]=tau[0][0]-tau[2][2];
	S[2]=tau[0][1];
	S[3]=2*tau[0][2];
	S[4]=tau[0][1]*U[2]/U[0]+tau[0][0]*U[1]/U[0]+tau[0][2]*U[3]/U[0]+NInitial::get_mu()/Pr*Hr;
}
void NMethods::stress(const int i, const int j, const double ***U, double (&tau)[3][3], char *bz, char *br,
	const double *hr, const double *hz)
{
	//char *bz="inout",*br="slip";
	double **gradUij, *U1, gradVelosity[2][3], rij;

	gradUij = new double *[2];
	for(int i = 0; i < 2; i++)
		gradUij[i] = new double [NInitial::getNcomp()];
	U1 = new double [NInitial::getNcomp()];
	NBoundaryCond *bound = new NBoundaryCond;

	get_grad1Visc(i,j,U,gradUij,bz,br,hr,hz);

	rij=rij1(i,hr);

	if((i>-1) && (i<NInitial::get_xmax()) && (j>-1) && (j<NInitial::get_ymax()))
	{
		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U[i][j][0]*(gradUij[k][l+1]-gradUij[k][0]*U[i][j][l+1]/U[i][j][0]);//k - компонента градиента (r,z)
	//l - компонента скорости (r,z,phi)

		tau[0][0]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[0][0]-
			gradVelosity[1][1]-U[i][j][1]/(U[i][j][0]*rij));
		tau[1][1]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[1][1]-
			gradVelosity[0][0]-U[i][j][1]/(U[i][j][0]*rij));
		tau[2][2]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*U[i][j][1]/(U[i][j][0]*rij)-
			gradVelosity[1][1]-gradVelosity[0][0]);
		tau[1][0]=tau[0][1]=NInitial::get_mu()*(gradVelosity[0][1]+gradVelosity[1][0]);
		tau[1][2]=NInitial::get_mu()*(gradVelosity[1][2]);
		tau[0][2]=NInitial::get_mu()*(gradVelosity[0][2]-U[i][j][3]/(U[i][j][0]*rij));
	}
	else if(i==-1)
	{
		if(br=="slip")
			bound->slipinnerVisc(U[0][j],U1);
		else
			{/*r-conditions*/}

		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента скорости (r,z)
	
		tau[0][0]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[0][0]-gradVelosity[1][1]-
			U1[1]/(U1[0]*rij));
		tau[1][1]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[1][1]-
			gradVelosity[0][0]-U1[1]/(U1[0]*rij));
		tau[2][2]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*U1[1]/(U1[0]*rij)-
			gradVelosity[1][1]-gradVelosity[0][0]);
		tau[1][0]=tau[0][1]=NInitial::get_mu()*(gradVelosity[0][1]+gradVelosity[1][0]);
		tau[1][2]=NInitial::get_mu()*(gradVelosity[1][2]);
		tau[0][2]=NInitial::get_mu()*(gradVelosity[0][2]-U1[3]/(U1[0]*rij));
	}
	else if(i==NInitial::get_xmax())
	{
		if(br=="slip")
			bound->slipVisc(U[NInitial::get_xmax()-1][j],U1,"r");
		else
			{/*r-conditions*/}

		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента скорости (r,z)
	
		tau[0][0]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[0][0]-
			gradVelosity[1][1]-U1[1]/(U1[0]*rij));
		tau[1][1]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[1][1]-
			gradVelosity[0][0]-U1[1]/(U1[0]*rij));
		tau[2][2]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*U1[1]/(U1[0]*rij)-
			gradVelosity[1][1]-gradVelosity[0][0]);
		tau[1][0]=tau[0][1]=NInitial::get_mu()*(gradVelosity[0][1]+gradVelosity[1][0]);
		tau[1][2]=NInitial::get_mu()*(gradVelosity[1][2]);
		tau[0][2]=NInitial::get_mu()*(gradVelosity[0][2]-U1[3]/(U1[0]*rij));
	}
	else if(j==-1)
	{
		if(bz=="slip")
			bound->slipVisc(U[i][0],U1,"z");
		else if(bz=="per")
			bound->periodic(i,0,const_cast<double***>(U),U1,"z");
		else 
			bound->outflow(U[i][0],U1);

		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента скорости (r,z)
	
		tau[0][0]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[0][0]-
			gradVelosity[1][1]-U1[1]/(U1[0]*rij));
		tau[1][1]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[1][1]-
			gradVelosity[0][0]-U1[1]/(U1[0]*rij));
		tau[2][2]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*U1[1]/(U1[0]*rij)-
			gradVelosity[1][1]-gradVelosity[0][0]);
		tau[1][0]=tau[0][1]=NInitial::get_mu()*(gradVelosity[0][1]+gradVelosity[1][0]);
		tau[1][2]=NInitial::get_mu()*(gradVelosity[1][2]);
		tau[0][2]=NInitial::get_mu()*(gradVelosity[0][2]-U1[3]/(U1[0]*rij));
	}
	else if(j==NInitial::get_ymax())
	{
		if(bz=="slip")
			bound->slipVisc(U[i][NInitial::get_ymax()-1],U1,"z");
		else if(bz=="per")
			bound->periodic(i,NInitial::get_ymax()-1,const_cast<double***>(U),U1,"z");
		else 
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
			NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);

		for(int k=0;k<2;k++)
			for(int l=0;l<3;l++)
				gradVelosity[k][l]=1/U1[0]*(gradUij[k][l+1]-gradUij[k][0]*U1[l+1]/U1[0]);//k - компонента скорости (r,z)
	
		tau[0][0]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[0][0]-
			gradVelosity[1][1]-U1[1]/(U1[0]*rij));
		tau[1][1]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*gradVelosity[1][1]-
			gradVelosity[0][0]-U1[1]/(U1[0]*rij));
		tau[2][2]=0.66666666666666666666666666666667*NInitial::get_mu()*(2*U1[1]/(U1[0]*rij)-
			gradVelosity[1][1]-gradVelosity[0][0]);
		tau[1][0]=tau[0][1]=NInitial::get_mu()*(gradVelosity[0][1]+gradVelosity[1][0]);
		tau[1][2]=NInitial::get_mu()*(gradVelosity[1][2]);
		tau[0][2]=NInitial::get_mu()*(gradVelosity[0][2]-U1[3]/(U1[0]*rij));
	}


	for(int i = 0; i < 2; i++)
		delete []gradUij[i];
	delete []gradUij;
	delete []U1;
	delete bound;
}
void NMethods::Fzvisc(const double *U, double *F, double tau[3][3], const double Hz)
{		
	double Pr;

	Pr=NInitial::getCp()*NInitial::get_mu()/NInitial::get_kappa();
	F[0]=0;
	F[1]=tau[1][0];
	F[2]=tau[1][1];
	F[3]=tau[1][2];
	F[4]=tau[1][1]*U[2]/U[0]+tau[1][0]*U[1]/U[0]+tau[1][2]*U[3]/U[0]+NInitial::get_mu()/Pr*Hz;
}
void NMethods::Frvisc(const double *U, double *F, double tau[3][3], const double Hr)
{		
	double Pr;

	Pr=NInitial::getCp()*NInitial::get_mu()/NInitial::get_kappa();
	F[0]=0;
	F[1]=tau[0][0];
	F[2]=tau[0][1];
	F[3]=tau[0][2];
	F[4]=tau[0][1]*U[2]/U[0]+tau[0][0]*U[1]/U[0]+tau[0][2]*U[3]/U[0]+NInitial::get_mu()/Pr*Hr;
}
void NMethods::limi(const int i, const int j, const double ***U, double *psi, char *dir, char *bz, char *br,
	const double *hr, const double *hz, const char *dirn)
{
	double **grUi, *delta2, *maxUj, *minUj, *Uprom, *Uadj, *Umax, *Umin;
	
	grUi = new double *[2];
	for(int i = 0; i < 2; i++)
		grUi[i] = new double[NInitial::getNcomp()];

	delta2 = new double[NInitial::getNcomp()];
	maxUj = new double[NInitial::getNcomp()];
	minUj = new double[NInitial::getNcomp()];
	Umax = new double[NInitial::getNcomp()];
	Umin = new double[NInitial::getNcomp()];

	get_grad1(i,j,const_cast<double***>(U),grUi,bz,br,hr,hz);

	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		if(dir=="r")
			if(dirn=="+")
				delta2[k]=grUi[0][k]*hr[i]*.5;
			else
				delta2[k]=grUi[0][k]*hr[i]*.5;
		else if(dir=="z")
			if(dirn=="+")
				delta2[k]=grUi[1][k]*hz[j]*.5;
			else
				delta2[k]=grUi[1][k]*hz[j]*.5;
		else
		{
			cout<<"err in lim dir";
			exit(1);
		}
	}

	vector<int> sizes(2);
	sizes[0] = NInitial::get_xmax();
	sizes[1] = NInitial::get_ymax();

	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		maxUj[k]=max(U[i-1+DUMMY_NUM2][j+DUMMY_NUM2][k], U[i+1+DUMMY_NUM2][j+DUMMY_NUM2][k]);
		minUj[k]=min(U[i-1+DUMMY_NUM2][j+DUMMY_NUM2][k], U[i+1+DUMMY_NUM2][j+DUMMY_NUM2][k]);
	}


	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		maxUj[k]=max(maxUj[k], U[i+DUMMY_NUM2][j+1+DUMMY_NUM2][k]);
		minUj[k]=min(minUj[k], U[i+DUMMY_NUM2][j+1+DUMMY_NUM2][k]);
	}

	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		maxUj[k]=max(maxUj[k], U[i+DUMMY_NUM2][j-1+DUMMY_NUM2][k]);
		minUj[k]=min(minUj[k], U[i+DUMMY_NUM2][j-1+DUMMY_NUM2][k]);
	}

	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		Umax[k]=max(U[i+DUMMY_NUM2][j+DUMMY_NUM2][k], maxUj[k]);
		Umin[k]=min(U[i+DUMMY_NUM2][j+DUMMY_NUM2][k], minUj[k]);
	}

	const double eps = .1e-5, omega = .1e-4;

	for(int k=0; k <= NInitial::getNcomp() - 1; k++)
	{
		if(abs(delta2[k]) < eps)
			delta2[k] = sign(delta2[k]) * (abs(delta2[k]) + omega);
			
		if(delta2[k] > 0)
			psi[k]=min(double(1), (Umax[k] - U[i+DUMMY_NUM2][j+DUMMY_NUM2][k])/delta2[k]);
		else if(delta2[k] < 0)
			psi[k]=min(double(1), (Umin[k] - U[i+DUMMY_NUM2][j+DUMMY_NUM2][k])/delta2[k]);
		else
			psi[k]=1;
	}

	for(int i = 0; i < 2; i++)
		delete []grUi[i];
	delete []grUi;

	delete []delta2;
	delete []maxUj;
	delete []minUj;
	delete []Umax;
	delete []Umin;
}
void NMethods::limibound(const int i, const int j, const double ***U, double *psi, char *dir, char *bz, char *br,
		const double *hr, const double *hz, const char *dirn)
{
	double *Uprom, *U1, *U2, *U3, *Umaxj, *Uminj, *Umax, *Umin,
		**grUi, *delta2, *U4;

	grUi = new double *[2];
	for(int i = 0; i < 2; i++)
		grUi[i] = new double[NInitial::getNcomp()];

	delta2 = new double[NInitial::getNcomp()];
	U1 = new double[NInitial::getNcomp()];
	U2 = new double[NInitial::getNcomp()];
	Uprom = new double[NInitial::getNcomp()];
	U3 = new double[NInitial::getNcomp()];
	Umax = new double[NInitial::getNcomp()];
	Umin = new double[NInitial::getNcomp()];
	Umaxj = new double[NInitial::getNcomp()];
	Uminj = new double[NInitial::getNcomp()];
	U4 = new double[NInitial::getNcomp()];

	NMgdBoundaryCond *bound = new NMgdBoundaryCond;

	if((i==-1) && (j==0))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+1][l];

			bound->slipinner(Uprom,U1);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+2][j][l];

			bound->slipinner(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j][l];

			bound->slipinner(Uprom,U3);

		}
		else
		{
			bound->periodic(i+1,j+1,const_cast<double***>(U),U1,"r");
			bound->periodic(i+2,j,const_cast<double***>(U),U2,"r");
			bound->periodic(i+1,j,const_cast<double***>(U),U3,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i+1][j][k],U1[k]);
			Uminj[k]=min(U[i+1][j][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((i==-1) && (j==NInitial::get_ymax()-1))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+2][j][l];

			bound->slipinner(Uprom,U1);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-1][l];

			bound->slipinner(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j][l];

			bound->slipinner(Uprom,U3);

		}
		else
		{
			bound->periodic(i+2,j,const_cast<double***>(U),U1,"r");
			bound->periodic(i+1,j-1,const_cast<double***>(U),U2,"r");
			bound->periodic(i+1,j,const_cast<double***>(U),U3,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i+1][j][k],U1[k]);
			Uminj[k]=min(U[i+1][j][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==-1) && (i==0))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+1][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->slip(Uprom,U2,"z");

		}
		else if(bz=="per")
		{
			bound->periodic(i+1,j+1,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j+2,const_cast<double***>(U),U2,"z");
		}
		else
		{			
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+1][l];

			bound->outflow(Uprom,U1);
			
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][0][l];

			bound->outflow(Uprom,U2);

		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+1][l];

			bound->slipinner(Uprom,U3);

		}
		else
		{
			bound->periodic(i,j+1,const_cast<double***>(U),U3,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i][j+1][k],U1[k]);
			Uminj[k]=min(U[i][j+1][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==-1) && (i==NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->slip(Uprom,U2,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i-1,j+1,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j+2,const_cast<double***>(U),U2,"z");
		}
		else
		{			
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->outflow(Uprom,U1);
			
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][0][l];

			bound->outflow(Uprom,U2);

		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+1][l];

			bound->slip(Uprom,U3,"r");
		}
		else
			bound->periodic(i,j+1,const_cast<double***>(U),U3,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i][j+1][k],U1[k]);
			Uminj[k]=min(U[i][j+1][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==0) && (i==NInitial::get_xmax()))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->slip(Uprom,U1,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-2][j][l];

			bound->slip(Uprom,U2,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j][l];

			bound->slip(Uprom,U3,"r");

		}
		else
		{
			bound->periodic(i-1,j+1,const_cast<double***>(U),U1,"r");
			bound->periodic(i-2,j,const_cast<double***>(U),U2,"r");
			bound->periodic(i-1,j,const_cast<double***>(U),U3,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i-1][j][k],U1[k]);
			Uminj[k]=min(U[i-1][j][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()-1) && (i==NInitial::get_xmax()))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-1][l];

			bound->slip(Uprom,U1,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-2][j][l];

			bound->slip(Uprom,U2,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j][l];

			bound->slip(Uprom,U3,"r");

		}
		else
		{
			bound->periodic(i-1,j-1,const_cast<double***>(U),U1,"r");
			bound->periodic(i-2,j,const_cast<double***>(U),U2,"r");
			bound->periodic(i-1,j,const_cast<double***>(U),U3,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i-1][j][k],U1[k]);
			Uminj[k]=min(U[i-1][j][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()) && (i==NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-1][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-2][l];

			bound->slip(Uprom,U2,"z");

		}
		else if(bz=="per")
		{
			bound->periodic(i-1,j-1,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j-2,const_cast<double***>(U),U2,"z");
		}
		else
		{
			bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U2,hr);
		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-1][l];

			bound->slip(Uprom,U3,"r");

		}
		else
			bound->periodic(i,j-1,const_cast<double***>(U),U3,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i][j-1][k],U1[k]);
			Uminj[k]=min(U[i][j-1][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()) && (i==0))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-1][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-2][l];

			bound->slip(Uprom,U2,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i+1,j-1,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j-2,const_cast<double***>(U),U2,"z");
		}
		else
		{
			bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U2,hr);
		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-1][l];

			bound->slipinner(Uprom,U3);
		}
		else
			bound->periodic(i,j-1,const_cast<double***>(U),U3,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i][j-1][k],U1[k]);
			Uminj[k]=min(U[i][j-1][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((i==-1) && (j>0) && (j<NInitial::get_ymax()-1))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-1][l];

			bound->slipinner(Uprom,U1);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+2][j][l];

			bound->slipinner(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+1][l];

			bound->slipinner(Uprom,U3);

		}
		else
		{
			bound->periodic(i+1,j-1,const_cast<double***>(U),U1,"r");
			bound->periodic(i+2,j,const_cast<double***>(U),U2,"r");
			bound->periodic(i+1,j+1,const_cast<double***>(U),U3,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i+1][j][k],U1[k]);
			Uminj[k]=min(U[i+1][j][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((i==-1) && (j==-1))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+1][l];

			bound->slipinner(Uprom,U1);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+2][j+1][l];

			bound->slipinner(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+2][l];

			bound->slipinner(Uprom,U3);

		}
		else
		{
			bound->periodic(i+1,j+1,const_cast<double***>(U),U1,"r");
			bound->periodic(i+2,j+1,const_cast<double***>(U),U2,"r");
			bound->periodic(i+1,j+2,const_cast<double***>(U),U3,"r");
		}

		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+1][l];

			bound->slip(Uprom,U4,"z");
		}
		else if(bz=="per")
			bound->periodic(i+1,j+1,const_cast<double***>(U),U4,"z");
		else
		{	
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][0][l];

			bound->outflow(Uprom,U4);
		}


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}	
	else if((i==-1) && (j==NInitial::get_ymax()))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-1][l];

			bound->slipinner(Uprom,U1);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+2][j-1][l];

			bound->slipinner(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-2][l];

			bound->slipinner(Uprom,U3);
		}
		else
		{
			bound->periodic(i+1,j-1,const_cast<double***>(U),U1,"r");
			bound->periodic(i+2,j-1,const_cast<double***>(U),U2,"r");
			bound->periodic(i+1,j-2,const_cast<double***>(U),U3,"r");
		}

		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-1][l];

			bound->slip(Uprom,U4,"z");
		}
		else if(bz=="per")
			bound->periodic(i+1,j-1,const_cast<double***>(U),U4,"z");
		else
			bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
			NInitial::get_uzin(),NInitial::get_uphiin(),U4,hr);


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}			
	else if((j==-1) && (i>0) && (i<NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->slip(Uprom,U2,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+1][l];

			bound->slip(Uprom,U3,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i-1,j+1,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j+2,const_cast<double***>(U),U2,"z");
			bound->periodic(i+1,j+1,const_cast<double***>(U),U3,"z");
		}
		else
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->outflow(Uprom,U1);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][0][l];

			bound->outflow(Uprom,U2);
			bound->outflow(Uprom,U3);
		}


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i][j+1][k],U1[k]);
			Uminj[k]=min(U[i][j+1][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==-1) && (i==NInitial::get_xmax()))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->slip(Uprom,U1,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-2][j+1][l];

			bound->slip(Uprom,U2,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+2][l];

			bound->slip(Uprom,U3,"r");

		}
		else
		{
			bound->periodic(i-1,j+1,const_cast<double***>(U),U1,"r");
			bound->periodic(i-2,j+1,const_cast<double***>(U),U2,"r");
			bound->periodic(i-1,j+2,const_cast<double***>(U),U3,"r");
		}

		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->slip(Uprom,U4,"z");
		}
		else if(bz=="per")
			bound->periodic(i-1,j+1,const_cast<double***>(U),U4,"z");
		else
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->outflow(Uprom,U4);
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((i==NInitial::get_xmax()) && (j>0) && (j<NInitial::get_ymax()-1))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-1][l];

			bound->slip(Uprom,U1,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-2][j][l];

			bound->slip(Uprom,U2,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+1][l];

			bound->slip(Uprom,U3,"r");
		}
		else
		{
			bound->periodic(i-1,j-1,const_cast<double***>(U),U1,"r");
			bound->periodic(i-2,j,const_cast<double***>(U),U2,"r");
			bound->periodic(i-1,j+1,const_cast<double***>(U),U3,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i-1][j][k],U1[k]);
			Uminj[k]=min(U[i-1][j][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()) && (i>0) && (i<NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-1][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-2][l];

			bound->slip(Uprom,U2,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-1][l];

			bound->slip(Uprom,U3,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i-1,j-1,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j-2,const_cast<double***>(U),U2,"z");
			bound->periodic(i+1,j-1,const_cast<double***>(U),U3,"z");
		}
		else
		{
			bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U2,hr);
			bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U3,hr);
		}


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U[i][j-1][k],U1[k]);
			Uminj[k]=min(U[i][j-1][k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()) && (i==NInitial::get_xmax()))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-1][l];

			bound->slip(Uprom,U1,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-2][j-1][l];

			bound->slip(Uprom,U2,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-2][l];

			bound->slip(Uprom,U3,"r");
		}
		else
		{
			bound->periodic(i-1,j-1,const_cast<double***>(U),U1,"r");
			bound->periodic(i-2,j-1,const_cast<double***>(U),U2,"r");
			bound->periodic(i-1,j-2,const_cast<double***>(U),U3,"r");
		}

		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-1][l];

			bound->slip(Uprom,U4,"z");
		}
		else if(bz=="per")
			bound->periodic(i-1,j-1,const_cast<double***>(U),U4,"z");
		else
			bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
			NInitial::get_uzin(),NInitial::get_uphiin(),U4,hr);


		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((i==-2) && (j>0) && (j<NInitial::get_ymax()-1))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+3][j+1][l];

			bound->slipinner(Uprom,U1);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+4][j][l];

			bound->slipinner(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+3][j-1][l];

			bound->slipinner(Uprom,U3);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+2][j][l];

			bound->slipinner(Uprom,U4);
		}
		else
		{
			bound->periodic(i+3,j+1,const_cast<double***>(U),U1,"r");
			bound->periodic(i+4,j,const_cast<double***>(U),U2,"r");
			bound->periodic(i+3,j-1,const_cast<double***>(U),U3,"r");
			bound->periodic(i+2,j,const_cast<double***>(U),U4,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((i==NInitial::get_xmax()+1) && (j>0) && (j<NInitial::get_ymax()-1))
	{
		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-3][j+1][l];

			bound->slip(Uprom,U1,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-4][j][l];

			bound->slip(Uprom,U2,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-3][j-1][l];

			bound->slip(Uprom,U3,"r");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-2][j][l];

			bound->slip(Uprom,U4,"r");
		}
		else
		{
			bound->periodic(i-3,j+1,const_cast<double***>(U),U1,"r");
			bound->periodic(i-4,j,const_cast<double***>(U),U2,"r");
			bound->periodic(i-3,j-1,const_cast<double***>(U),U3,"r");
			bound->periodic(i-2,j,const_cast<double***>(U),U4,"r");
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==-2) && (i>0) && (i<NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+4][l];

			bound->slip(Uprom,U2,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+3][l];

			bound->slip(Uprom,U3,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+3][l];

			bound->slip(Uprom,U4,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i,j+2,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j+4,const_cast<double***>(U),U2,"z");
			bound->periodic(i-1,j+3,const_cast<double***>(U),U3,"z");
			bound->periodic(i+1,j+3,const_cast<double***>(U),U4,"z");
		}
		else
		{		
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->outflow(Uprom,U1);
			bound->outflow(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][0][l];

			bound->outflow(Uprom,U3);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][0][l];

			bound->outflow(Uprom,U4);
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()+1) && (i>0) && (i<NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-2][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-4][l];

			bound->slip(Uprom,U2,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-3][l];

			bound->slip(Uprom,U3,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-3][l];

			bound->slip(Uprom,U4,"z");

		}
		else if(bz=="per")
		{
			bound->periodic(i,j-2,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j-4,const_cast<double***>(U),U2,"z");
			bound->periodic(i-1,j-3,const_cast<double***>(U),U3,"z");
			bound->periodic(i+1,j-3,const_cast<double***>(U),U4,"z");
		}
		else
		{
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U2,hr);
			bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U3,hr);
			bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U4,hr);
		}

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==-2) && (i==0))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+4][l];

			bound->slip(Uprom,U2,"z");


			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j+3][l];

			bound->slip(Uprom,U4,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i,j+2,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j+4,const_cast<double***>(U),U2,"z");
			bound->periodic(i+1,j+3,const_cast<double***>(U),U4,"z");
		}
		else
		{	
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->outflow(Uprom,U1);
			bound->outflow(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][0][l];

			bound->outflow(Uprom,U4);
		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+3][l];

			bound->slipinner(Uprom,U3);
		}
		else	
			bound->periodic(i,j+3,const_cast<double***>(U),U3,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()+1) && (i==0))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-2][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-4][l];

			bound->slip(Uprom,U2,"z");


			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i+1][j-3][l];

			bound->slip(Uprom,U4,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i,j-2,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j-4,const_cast<double***>(U),U2,"z");
			bound->periodic(i+1,j-3,const_cast<double***>(U),U4,"z");
		}
		else
		{
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U2,hr);
			bound->inflow(i+1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U4,hr);
		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-3][l];

			bound->slipinner(Uprom,U3);
		}
		else	
			bound->periodic(i,j-3,const_cast<double***>(U),U3,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==-2) && (i==NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+4][l];

			bound->slip(Uprom,U2,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j+3][l];

			bound->slip(Uprom,U3,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i,j+2,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j+4,const_cast<double***>(U),U2,"z");
			bound->periodic(i-1,j+3,const_cast<double***>(U),U3,"z");
		}
		else
		{	
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+2][l];

			bound->outflow(Uprom,U1);
			bound->outflow(Uprom,U2);

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][0][l];

			bound->outflow(Uprom,U3);
		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j+3][l];

			bound->slip(Uprom,U4,"r");
		}
		else
			bound->periodic(i,j+3,const_cast<double***>(U),U4,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}
	else if((j==NInitial::get_ymax()+1) && (i==NInitial::get_xmax()-1))
	{
		if(bz=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-2][l];

			bound->slip(Uprom,U1,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-4][l];

			bound->slip(Uprom,U2,"z");

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i-1][j-3][l];

			bound->slip(Uprom,U3,"z");
		}
		else if(bz=="per")
		{
			bound->periodic(i,j-2,const_cast<double***>(U),U1,"z");
			bound->periodic(i,j-4,const_cast<double***>(U),U2,"z");
			bound->periodic(i-1,j-3,const_cast<double***>(U),U3,"z");
		}
		else
		{
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U2,hr);
			bound->inflow(i-1,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
				NInitial::get_uzin(),NInitial::get_uphiin(),U3,hr);
		}

		if(br=="slip")
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Uprom[l]=U[i][j-3][l];

			bound->slip(Uprom,U4,"r");
		}
		else
			bound->periodic(i,j-3,const_cast<double***>(U),U4,"r");

		for(int k=0;k<=NInitial::getNcomp()-1;k++)
		{
			Umaxj[k]=max(U4[k],U1[k]);
			Uminj[k]=min(U4[k],U1[k]);

			Umaxj[k]=max(Umaxj[k],U2[k]);
			Uminj[k]=min(Uminj[k],U2[k]);

			Uminj[k]=min(U3[k],Uminj[k]);
			Umaxj[k]=max(U3[k],Umaxj[k]);
		}
	}

	double *Ucell, *Ust;

	Ucell = new double[NInitial::getNcomp()];
	Ust = new double[NInitial::getNcomp()];

	if(i==-1)
	{
		if(j == -1)
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Ust[l]=U[i+1][j+1][l];

			double *UcellForCell = new double[NInitial::getNcomp()];//костыль))))))))

			if(bz=="slip")
				bound->slip(Ust,UcellForCell,"z");
			else
				bound->periodic(i+1,j+1,const_cast<double***>(U),UcellForCell,"z");

			double ***UcellForCellStrih;//ужасный костыль, за это надо увольн€ть сразу, но ничего не поделаешь))))))))
			
			UcellForCellStrih = new double**[NInitial::get_xmax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				UcellForCellStrih[i] = new double*[NInitial::get_ymax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					UcellForCellStrih[i][j] = new double[NInitial::getNcomp()];

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				UcellForCellStrih[0][0][l] = UcellForCell[l];

			if(br=="slip")
				bound->slipinner(UcellForCell,Ucell);
			else
				bound->periodic(i+1,j+1,const_cast<double***>(UcellForCellStrih),Ucell,"r");

			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					delete []UcellForCellStrih[i][j];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				delete []UcellForCellStrih[i];
			delete []UcellForCellStrih;
			delete []UcellForCell;
		}
		else if(j == NInitial::get_ymax())
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Ust[l]=U[i+1][j-1][l];

			double *UcellForCell = new double[NInitial::getNcomp()];//костыль))))))))

			if(bz=="slip")
				bound->slip(Ust,UcellForCell,"z");
			else
				bound->periodic(i+1,j-1,const_cast<double***>(U),UcellForCell,"z");

			double ***UcellForCellStrih;//ужасный костыль, за это надо увольн€ть сразу, но ничего не поделаешь))))))))
			
			UcellForCellStrih = new double**[NInitial::get_xmax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				UcellForCellStrih[i] = new double*[NInitial::get_ymax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					UcellForCellStrih[i][j] = new double[NInitial::getNcomp()];

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				UcellForCellStrih[0][j-1][l] = UcellForCell[l];

			if(br=="slip")
				bound->slipinner(UcellForCell,Ucell);
			else
				bound->periodic(i+1,j-1,const_cast<double***>(UcellForCellStrih),Ucell,"r");

			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					delete []UcellForCellStrih[i][j];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				delete []UcellForCellStrih[i];
			delete []UcellForCellStrih;
			delete []UcellForCell;
		}
		else
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Ust[l]=U[i+1][j][l];

			if(br=="slip")
				bound->slipinner(Ust,Ucell);
			else
				bound->periodic(i+1,j,const_cast<double***>(U),Ucell,"r");
		}
	}
	else if(i==NInitial::get_xmax())
	{
		if(j == -1)
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Ust[l]=U[i-1][j+1][l];

			double *UcellForCell = new double[NInitial::getNcomp()];//костыль))))))))

			if(bz=="slip")
				bound->slip(Ust,UcellForCell,"z");
			else
				bound->periodic(i-1,j+1,const_cast<double***>(U),UcellForCell,"z");

			double ***UcellForCellStrih;//ужасный костыль, за это надо увольн€ть сразу, но ничего не поделаешь))))))))
			
			UcellForCellStrih = new double**[NInitial::get_xmax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				UcellForCellStrih[i] = new double*[NInitial::get_ymax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					UcellForCellStrih[i][j] = new double[NInitial::getNcomp()];

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				UcellForCellStrih[i-1][j+1][l] = UcellForCell[l];

			if(br=="slip")
				bound->slipinner(UcellForCell,Ucell);
			else
				bound->periodic(i-1,j+1,const_cast<double***>(UcellForCellStrih),Ucell,"r");

			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					delete []UcellForCellStrih[i][j];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				delete []UcellForCellStrih[i];
			delete []UcellForCellStrih;
			delete []UcellForCell;
		}
		else if(j == NInitial::get_ymax())
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Ust[l]=U[i-1][j-1][l];

			double *UcellForCell = new double[NInitial::getNcomp()];//костыль))))))))

			if(bz=="slip")
				bound->slip(Ust,UcellForCell,"z");
			else
				bound->periodic(i-1,j-1,const_cast<double***>(U),UcellForCell,"z");

			double ***UcellForCellStrih;//ужасный костыль, за это надо увольн€ть сразу, но ничего не поделаешь))))))))
			
			UcellForCellStrih = new double**[NInitial::get_xmax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				UcellForCellStrih[i] = new double*[NInitial::get_ymax()];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					UcellForCellStrih[i][j] = new double[NInitial::getNcomp()];

			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				UcellForCellStrih[i-1][j-1][l] = UcellForCell[l];

			if(br=="slip")
				bound->slipinner(UcellForCell,Ucell);
			else
				bound->periodic(i-1,j-1,const_cast<double***>(UcellForCellStrih),Ucell,"r");

			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t j = 0; j < NInitial::get_ymax(); j++)
					delete []UcellForCellStrih[i][j];
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				delete []UcellForCellStrih[i];
			delete []UcellForCellStrih;
			delete []UcellForCell;
		}
		else
		{
			for(int l=0;l<=NInitial::getNcomp()-1;l++)
				Ust[l]=U[i-1][j][l];

			if(br=="slip")
				bound->slip(Ust,Ucell,"r");
			else
				bound->periodic(i-1,j,const_cast<double***>(U),Ucell,"r");
		}
	}
	else if(j==-1)
	{
		for(int l=0;l<=NInitial::getNcomp()-1;l++)
			Ust[l]=U[i][j+1][l];

		if(bz=="slip")
			bound->slip(Ust,Ucell,"z");
		else if(bz=="per")
			bound->periodic(i,j+1,const_cast<double***>(U),Ucell,"z");
		else
			bound->outflow(Ust,Ucell);
	}
	else if(j==NInitial::get_ymax())
	{
		for(int l=0;l<=NInitial::getNcomp()-1;l++)
			Ust[l]=U[i][j-1][l];

		if(bz=="slip")
			bound->slip(Ust,Ucell,"z");
		else if(bz=="per")
			bound->periodic(i,j-1,const_cast<double***>(U),Ucell,"z");
		else
			bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),NInitial::get_urin(),
			NInitial::get_uzin(),NInitial::get_uphiin(),Ucell,hr);
	}
	else
		for(int l=0;l<=NInitial::getNcomp()-1;l++)
			Ucell[l]=Ust[l];
		
	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		Umax[k]=max(Ucell[k],Umaxj[k]);
		Umin[k]=min(Ucell[k],Uminj[k]);
	}

	get_grad1(i,j,const_cast<double***>(U),grUi,bz,br,hr,hz);

	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		if(dir=="r")
		{	
			if(dirn=="+")
			{
				if(i>=0)
					delta2[k]=grUi[0][k]*hr[i]*.5;
				else if(i==-1)
					delta2[k]=grUi[0][k]*hr[0]*.5;
				else if(i==-2)
					delta2[k]=grUi[0][k]*hr[1]*.5;
				else if(i==-3)
					delta2[k]=grUi[0][k]*hr[2]*.5;
			}
			else
			{
				if(i>=0)
					delta2[k]=grUi[0][k]*hr[i]*.5;
				else if(i==-1)
					delta2[k]=grUi[0][k]*hr[0]*.5;
				else if(i==-2)
					delta2[k]=grUi[0][k]*hr[1]*.5;
				else if(i==-3)
					delta2[k]=grUi[0][k]*hr[2]*.5;
			}
		}
		else if(dir=="z")
		{	
			if(dirn=="+")
			{
				if(j>=0)
					delta2[k]=grUi[1][k]*hz[j]*.5;
				else if(j==-1)
					delta2[k]=grUi[1][k]*hz[0]*.5;
				else if(j==-2)
					delta2[k]=grUi[1][k]*hz[1]*.5;
				else if(j==-3)
					delta2[k]=grUi[1][k]*hz[2]*.5;
			}
			else
			{
				if(j>=0)
					delta2[k]=grUi[1][k]*hz[j]*.5;
				else if(j==-1)
					delta2[k]=grUi[1][k]*hz[0]*.5;
				else if(j==-2)
					delta2[k]=grUi[1][k]*hz[1]*.5;
				else if(j==-3)
					delta2[k]=grUi[1][k]*hz[2]*.5;
			}
		}
		else
		{
			cout<<"err in lim dir";
			exit(1);
		}
	}
		
	const double eps=.1e-5, omega=.1e-4;

	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		if(abs(delta2[k]) < eps)
			delta2[k]=sign(delta2[k])*(abs(delta2[k])+omega);

		if(delta2[k]>0)
			psi[k]=min(double(1),(Umax[k]-Ucell[k])/delta2[k]);
		else if(delta2[k]<0)
			psi[k]=min(double(1),(Umin[k]-Ucell[k])/delta2[k]);
		else
			psi[k]=1;
	}

	delete []Ucell;
	delete []Ust;

	for(int i = 0; i < 2; i++)
		delete []grUi[i];
	delete []grUi;

	delete []delta2;
	delete []U1;
	delete []U2;
	delete []Uprom;
	delete []U3;
	delete []Umax;
	delete []Umin;
	delete []Umaxj;
	delete []Uminj;
	delete []U4;

	delete bound;

}
void NMethods::limiter(const int i, const int j, double ***U, double *psi, char *bz, char *br,
		const double *hr, const double *hz, const char *dirn)
{
	double *psi2, *psi3, *psi4;

	psi2 = new double[NInitial::getNcomp()];
	psi3 = new double[NInitial::getNcomp()];
	psi4 = new double[NInitial::getNcomp()];

	limi(i+1,j,const_cast<const double***>(U),psi,"r",bz,br,hr,hz,dirn);//i=[0;xmax-1),j=[0;ymax-1]
	
	limi(i-1,j,const_cast<const double***>(U),psi2,"r",bz,br,hr,hz,dirn);//i=(0;xmax-1],j=[0;ymax-1]

	limi(i,j+1,const_cast<const double***>(U),psi3,"z",bz,br,hr,hz,dirn);//i=[0;xmax-1],j=[0;ymax-1)

	limi(i,j-1,const_cast<const double***>(U),psi4,"z",bz,br,hr,hz,dirn);//i=[0;xmax-1],j=(0;ymax-1]


	for(int k=0;k<=NInitial::getNcomp()-1;k++)
	{
		psi[k]=min(psi2[k],psi[k]);
		psi2[k]=min(psi2[k],psi3[k]);
		psi2[k]=min(psi[k],psi2[k]);
		psi[k]=min(psi2[k],psi4[k]);
	}

	delete []psi2;
	delete []psi3;
	delete []psi4;
}
double NMethods::rij1(const int &i, const double *hr)//h0 - шаг первой €чейки слева
{//Si = h[0] * (1-pow(RFAKTOR,i))/(1-RFAKTOR), i=[0,n]
	double sum;

	if((i >= 0) & (i < NInitial::get_xmax()))
	{
		sum = NInitial::getR1();
		for(int k = 0; k < i; k++)
			sum = sum + hr[k];
		sum = sum + .5 * hr[i];
	}
	else if(i >= NInitial::get_xmax())
	{
		//int n;
		//n = i - xmax + 1;
		sum = NInitial::getR1();
		for(int k = 0; k < 2 * NInitial::get_xmax() - i - 1; k++)
			sum = sum + hr[k];
		sum = sum + .5 * hr[2 * NInitial::get_xmax() - i - 1];
	}
	else
	{
		sum = NInitial::getR1();
		for(int k = -1; k > i; k--)
			sum = sum + hr[abs(k) - 1];
		sum = sum + .5 * hr[abs(i) - 1];
	}

	return sum;
}
double NMethods::sound(const double *U, const char *dir)
{
	double b2;
	b2 = 1/U[0] * (pow(U[5], 2) + pow(U[6], 2) + pow(U[7], 2));
	//greater magnetogydrodynamic wave
	if(!strcmp(dir, "r"))
	{
		return sqrt(.5 * (pow(sound(U), 2) + b2 + 
			sqrt(pow(pow(sound(U), 2) + b2, 2) - 4. * pow(sound(U), 2) * pow(U[5], 2)/U[0])));
	}
	else
	{
		return sqrt(.5 * (pow(sound(U), 2) + b2 + 
			sqrt(pow(pow(sound(U), 2) + b2, 2) - 4. * pow(sound(U), 2) * pow(U[6], 2)/U[0])));
	}
}



#endif //METHODS_CPP