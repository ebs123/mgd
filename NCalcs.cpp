#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\NCalcs.h"
#else
#include "includes/NCalcs.h"
#endif


NCalcs::NCalcs(void)
{
}
NCalcs::~NCalcs(void)
{
}
void NCalcs::resid(double ***U, double ***Re,	const double *hr, const double *hz, double ***Rv)
{
	double *U1, *Ul, *Ur, *Fr1, *Fr2, *Fz1, *Fz2, *S,
		dlr, dlz, **gradUr, **gradUl, *psiL, rij, *Frv1, *Frv2, *Fzv1, *Fzv2,
		*psiR, *Sv, gradH[2];
	char *bz = "per", *br = "slip";

	U1 = new double[NInitial::getNcomp()];
	Ul = new double[NInitial::getNcomp()];
	Ur = new double[NInitial::getNcomp()];
	Fr1 = new double[NInitial::getNcomp()];
	Fr2 = new double[NInitial::getNcomp()];
	Fz1 = new double[NInitial::getNcomp()];
	Fz2 = new double[NInitial::getNcomp()];
	S = new double[NInitial::getNcomp()];
	gradUr = new double*[2];
	for(size_t i = 0; i < 2; i++)
		gradUr[i] = new double[NInitial::getNcomp()];
	gradUl = new double*[2];
	for(size_t i = 0; i < 2; i++)
		gradUl[i] = new double[NInitial::getNcomp()];
	psiL = new double[NInitial::getNcomp()];
	Frv1 = new double[NInitial::getNcomp()];
	Frv2 = new double[NInitial::getNcomp()];
	Fzv1 = new double[NInitial::getNcomp()];
	Fzv2 = new double[NInitial::getNcomp()];
	psiR = new double[NInitial::getNcomp()];
	//Sv = new double[NInitial::getNcomp()];

	NMgdMethods *methods = new NMgdMethods;
	NMgdBoundaryCond *bound = new NMgdBoundaryCond;
	vector<int> sizes(2);
	sizes[0] = NInitial::get_xmax();
	sizes[1] = NInitial::get_ymax();
	NArrPacker<double> *U_pack = new NArrPacker<double>(DUMMY_NUM2, sizes, NInitial::getPeriodicName(), NInitial::getSlipName(), (void*)U);
	double*** U_pack_arr = (double***)U_pack->getPackArr();


	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
		{
			try
			{
			rij = methods->rij1(i, hr);
			//double tauR[3][3], tauL[3][3], tau[3][3], *Uaver, gradHL[2], gradHR[2];

			//Uaver = new double[NInitial::getNcomp()];

			if((i != 0) & (j != 0) & (j != NInitial::get_ymax() - 1) & (i != NInitial::get_xmax() - 1))
			{
				methods->get_grad1(i, j, U_pack_arr, gradUr, bz, br, hr, hz);
				methods->get_grad1(i - 1, j, U_pack_arr, gradUl, bz, br, hr, hz);

				methods->limiter(i-1,j,U_pack_arr,psiL,bz,br,hr,hz,"-");
				methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");

				//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
				//methods->stress(i-1,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

				//for(int l=0;l<3;l++)
				//	for(int m=0;m<3;m++)
				//		tau[l][m]=(tauL[l][m]*hr[i]+tauR[l][m]*hr[i-1])/(hr[i-1]+hr[i]);

				for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				{
					Ul[k]=U[i-1][j][k]+psiL[k]*gradUl[0][k]*hr[i-1]*.5;
					Ur[k]=U[i][j][k]-psiR[k]*gradUr[0][k]*hr[i]*.5;
				}
	
				//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
				//methods->getGradH(i-1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

				//for(int k=0;k<2;k++)
				//	gradH[k]=(gradHL[k]*hr[i]+gradHR[k]*hr[i-1])/(hr[i-1]+hr[i]);
				
				//methods->Frvisc(Uaver,Frv1,tau,gradH[0]);
				methods->NMgdMethods::Roe(Ur,Ul,Fr1,"r","-");

				if(isnan(Fr1))
				{
					NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr1 is nan!");
					exit(1);
				}

				methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
				methods->get_grad1(i+1,j,U_pack_arr,gradUr,bz,br,hr,hz);

				methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
				methods->limiter(i+1,j,U_pack_arr,psiR,bz,br,hr,hz,"+");

				//methods->stress(i+1,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
				//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

				//for(int l=0;l<3;l++)
				//	for(int m=0;m<3;m++)
				//		tau[l][m]=(tauL[l][m]*hr[i+1]+tauR[l][m]*hr[i])/(hr[i+1]+hr[i]);

				for(int k=0; k <= NInitial::getNcomp() - 1; k++)
				{
					Ul[k]=U[i][j][k]+psiL[k]*gradUl[0][k]*hr[i]*.5;
					Ur[k]=U[i+1][j][k]-psiR[k]*gradUr[0][k]*hr[i+1]*.5;

					//Uaver[k]=(U[i+1][j][k]*hr[i]+U[i][j][k]*hr[i+1])/(hr[i+1]+hr[i]);
				}

				//methods->getGradH(i+1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
				//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

				//for(int k=0;k<2;k++)
				//	gradH[k]=(gradHL[k]*hr[i+1]+gradHR[k]*hr[i])/(hr[i+1]+hr[i]);

				//methods->Frvisc(Uaver,Frv2,tau,gradH[0]);
				methods->NMgdMethods::Roe(Ur,Ul,Fr2,"r","+");

				if(isnan(Fr2))
				{
					NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr2 is nan!");
					exit(1);
				}

				methods->get_grad1(i,j,U_pack_arr,gradUr,bz,br,hr,hz);
				methods->get_grad1(i,j-1,U_pack_arr,gradUl,bz,br,hr,hz);

				methods->limiter(i,j-1,U_pack_arr,psiL,bz,br,hr,hz,"-");
				methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");

				//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
				//methods->stress(i,j-1,const_cast<const double***>(U),tauL,bz,br,hr,hz);

				//for(int l=0;l<3;l++)
				//	for(int m=0;m<3;m++)
				//		tau[l][m]=(tauL[l][m]*hz[j]+tauR[l][m]*hz[j-1])/(hz[j-1]+hz[j]);

				for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				{
					Ul[k]=U[i][j-1][k]+psiL[k]*gradUl[1][k]*hz[j-1]*.5;
					Ur[k]=U[i][j][k]-psiR[k]*gradUr[1][k]*hz[j]*.5;

					//Uaver[k]=(U[i][j-1][k]*hz[j]+U[i][j][k]*hz[j-1])/(hz[j-1]+hz[j]);
				}

				//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
				//methods->getGradH(i,j-1,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

				//for(int k=0;k<2;k++)
				//	gradH[k]=(gradHL[k]*hz[j]+gradHR[k]*hz[j-1])/(hz[j-1]+hz[j]);

				//methods->Fzvisc(Uaver,Fzv1,tau,gradH[1]);
				methods->NMgdMethods::Roe(Ur,Ul,Fz1,"z","-");

				if(isnan(Fz1))
				{
					NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz1 is nan!");
					exit(1);
				}

				methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
				methods->get_grad1(i,j+1,U_pack_arr,gradUr,bz,br,hr,hz);

				methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
				methods->limiter(i,j+1,U_pack_arr,psiR,bz,br,hr,hz,"+");

				//methods->stress(i,j+1,const_cast<const double***>(U),tauR,bz,br,hr,hz);
				//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

				//for(int l=0;l<3;l++)
				//	for(int m=0;m<3;m++)
				//		tau[l][m]=(tauL[l][m]*hz[j+1]+tauR[l][m]*hz[j])/(hz[j+1]+hz[j]);

				for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				{
					Ul[k]=U[i][j][k]+psiL[k]*gradUl[1][k]*hz[j]*.5;
					Ur[k]=U[i][j+1][k]-psiR[k]*gradUr[1][k]*hz[j+1]*.5;

					//Uaver[k]=(U[i][j+1][k]*hz[j]+U[i][j][k]*hz[j+1])/(hz[j+1]+hz[j]);
				}

				//methods->getGradH(i,j+1,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
				//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

				//for(int k=0;k<2;k++)
				//	gradH[k]=(gradHL[k]*hz[j+1]+gradHR[k]*hz[j])/(hz[j+1]+hz[j]);

				//methods->Fzvisc(Uaver,Fzv2,tau,gradH[1]);
				methods->NMgdMethods::Roe(Ur,Ul,Fz2,"z","+");

				if(isnan(Fz2))
				{
					NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz2 is nan!");
					exit(1);
				}

			}
			else if((i == 0) || (i == NInitial::get_xmax() - 1))
			{
				double *Ust;

				Ust = new double[NInitial::getNcomp()];

				for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					Ust[k] = U[i][j][k];

				if(br == "slip")
					if(i == 0)
						bound->slipinner(Ust, U1);
					else
						bound->slip(Ust, U1, "r");
				else
					bound->periodic(i,j,U,U1,"r");

				if(i==0)
				{
					methods->get_grad1(i-1,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i-1,j,U_pack_arr,psiL,bz,br,hr,hz,"-");
					methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");

					//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i-1,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);
			
					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]+tauR[l][m])*.5;

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U1[k]+psiL[k]*gradUl[0][k]*hr[i]*.5;
						Ur[k]=U[i][j][k]-psiR[k]*gradUr[0][k]*hr[i]*.5;

						//Uaver[k]=(U1[k]+U[i][j][k])*.5;
					}
			
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i-1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]+gradHR[k])*.5;

					//methods->Frvisc(Uaver,Frv1,tau,gradH[0]);
					methods->NMgdMethods::Roe(Ur,Ul,Fr1,"r","-");
			
					if(isnan(Fr1))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr1 is nan!");
						exit(1);
					}

					methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i+1,j,U_pack_arr,gradUr,bz,br,hr,hz);
			
					methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
					methods->limiter(i+1,j,U_pack_arr,psiR,bz,br,hr,hz,"+");
			
					//methods->stress(i+1,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hr[i+1]+tauR[l][m]*hr[i])/(hr[i+1]+hr[i]);
			
					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j][k]+psiL[k]*gradUl[0][k]*hr[i]*.5;
						Ur[k]=U[i+1][j][k]-psiR[k]*gradUr[0][k]*hr[i+1]*.5;

						//Uaver[k]=(U[i+1][j][k]*hr[i]+U[i][j][k]*hr[i+1])/(hr[i+1]+hr[i]);
					}

					//methods->getGradH(i+1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]*hr[i+1]+gradHR[k]*hr[i])/(hr[i+1]+hr[i]);

					//methods->Frvisc(Uaver,Frv2,tau,gradH[0]);
					methods->NMgdMethods::Roe(Ur,Ul,Fr2,"r","+");

					if(isnan(Fr2))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr2 is nan!");
						exit(1);
					}

					//double *traceData = new double[NInitial::getNcomp() + 2];
					//traceData[NInitial::getNcomp()] = i;
					//traceData[NInitial::getNcomp() + 1] = j;
					//for(int l = 0; l < NInitial::getNcomp(); l++)
					//	traceData[l] = Ul[l];
					//NTracer::traceToFile(traceData, "double", NInitial::getNcomp() + 2);
					//for(int l = 0; l < NInitial::getNcomp(); l++)
					//	traceData[l] = Ur[l];
					//NTracer::traceToFile(traceData, "double", NInitial::getNcomp() + 2);
					//delete []traceData;
				}
				else
				{
					methods->get_grad1(i-1,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i-1,j,U_pack_arr,psiL,bz,br,hr,hz,"-");
					methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");

					//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i-1,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hr[i]+tauR[l][m]*hr[i-1])/(hr[i-1]+hr[i]);

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i-1][j][k]+psiL[k]*gradUl[0][k]*hr[i-1]*.5;
						Ur[k]=U[i][j][k]-psiR[k]*gradUr[0][k]*hr[i]*.5;

						//Uaver[k]=(U[i-1][j][k]*hr[i]+U[i][j][k]*hr[i-1])/(hr[i-1]+hr[i]);
					}

					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i-1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]*hr[i]+gradHR[k]*hr[i-1])/(hr[i-1]+hr[i]);

					//methods->Frvisc(Uaver,Frv1,tau,gradH[0]);
					methods->NMgdMethods::Roe(Ur,Ul,Fr1,"r","-");

					if(isnan(Fr1))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr1 is nan!");
						exit(1);
					}

					methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i+1,j,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
					methods->limiter(i+1,j,U_pack_arr,psiR,bz,br,hr,hz,"+");

					//methods->stress(i+1,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]+tauR[l][m])*.5;

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j][k]+psiL[k]*gradUl[0][k]*hr[i]*.5;
						Ur[k]=U1[k]-psiR[k]*gradUr[0][k]*hr[i]*.5;

						//Uaver[k]=(U1[k]+U[i][j][k])*.5;
					}

					//methods->getGradH(i+1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]+gradHR[k])*.5;

					//methods->Frvisc(Uaver,Frv2,tau,gradH[0]);
					methods->NMgdMethods::Roe(Ur,Ul,Fr2,"r","+");

					if(isnan(Fr2))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr2 is nan!");
						exit(1);
					}

				}

				if((j != 0) & (j != NInitial::get_ymax() - 1))
				{
					methods->get_grad1(i,j-1,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j-1,U_pack_arr,psiL,bz,br,hr,hz,"-");
					methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");

					//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j-1,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hz[j]+tauR[l][m]*hz[j-1])/(hz[j-1]+hz[j]);

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j-1][k]+psiL[k]*gradUl[1][k]*hz[j-1]*.5;
						Ur[k]=U[i][j][k]-psiR[k]*gradUr[1][k]*hz[j]*.5;

						//Uaver[k]=(U[i][j-1][k]*hz[j]+U[i][j][k]*hz[j-1])/(hz[j-1]+hz[j]);
					}

					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j-1,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]*hz[j]+gradHR[k]*hz[j-1])/(hz[j-1]+hz[j]);

					//methods->Fzvisc(Uaver,Fzv1,tau,gradH[1]);
					methods->NMgdMethods::Roe(Ur,Ul,Fz1,"z","-");

					if(isnan(Fz1))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz1 is nan!");
						exit(1);
					}

					methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j+1,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
					methods->limiter(i,j+1,U_pack_arr,psiR,bz,br,hr,hz,"+");

					//methods->stress(i,j+1,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hz[j+1]+tauR[l][m]*hz[j])/(hz[j+1]+hz[j]);

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j][k]+psiL[k]*gradUl[1][k]*hz[j]*.5;
						Ur[k]=U[i][j+1][k]-psiR[k]*gradUr[1][k]*hz[j+1]*.5;

						//Uaver[k]=(U[i][j+1][k]*hz[j]+U[i][j][k]*hz[j+1])/(hz[j+1]+hz[j]);
					}

					//methods->getGradH(i,j+1,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]*hz[j+1]+gradHR[k]*hz[j])/(hz[j+1]+hz[j]);

					//methods->Fzvisc(Uaver,Fzv2,tau,gradH[1]);
					methods->NMgdMethods::Roe(Ur,Ul,Fz2,"z","+");

					if(isnan(Fz2))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz2 is nan!");
						exit(1);
					}

				}

			delete []Ust;
			}

			if((j == 0) | (j == NInitial::get_ymax() - 1))
			{
				double *Ust;

				Ust = new double[NInitial::getNcomp()];

				for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					Ust[k]=U[i][j][k];

				if(bz == "slip")
					bound->slip(Ust,U1,"z");
				else if(bz == "per")
					bound->periodic(i,j,U,U1,"z");
				else
				{
					if(j == 0)
						bound->outflow(Ust,U1);
					else
						bound->inflow(i,NInitial::get_pin(),NInitial::get_roin(),
						NInitial::get_urin(),NInitial::get_uzin(),NInitial::get_uphiin(),U1,hr);
				}

				if((i != 0) & (i != NInitial::get_xmax() - 1))
				{
					methods->get_grad1(i-1,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j,U_pack_arr,gradUr,bz,br,hr,hz);
			
					methods->limiter(i-1,j,U_pack_arr,psiL,bz,br,hr,hz,"-");
					methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");
		
					//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i-1,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hr[i]+tauR[l][m]*hr[i-1])/(hr[i-1]+hr[i]);

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i-1][j][k]+psiL[k]*gradUl[0][k]*hr[i-1]*.5;
						Ur[k]=U[i][j][k]-psiR[k]*gradUr[0][k]*hr[i]*.5;

						//Uaver[k]=(U[i-1][j][k]*hr[i]+U[i][j][k]*hr[i-1])/(hr[i-1]+hr[i]);
					}

					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i-1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]*hr[i]+gradHR[k]*hr[i-1])/(hr[i-1]+hr[i]);

					//methods->Frvisc(Uaver,Frv1,tau,gradH[0]);
					methods->NMgdMethods::Roe(Ur,Ul,Fr1,"r","-");

					if(isnan(Fr1))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr1 is nan!");
						exit(1);
					}

					methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i+1,j,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
					methods->limiter(i+1,j,U_pack_arr,psiR,bz,br,hr,hz,"+");

					//methods->stress(i+1,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hr[i+1]+tauR[l][m]*hr[i])/(hr[i+1]+hr[i]);

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j][k]+psiL[k]*gradUl[0][k]*hr[i]*.5;
						Ur[k]=U[i+1][j][k]-psiR[k]*gradUr[0][k]*hr[i+1]*.5;

						//Uaver[k]=(U[i+1][j][k]*hr[i]+U[i][j][k]*hr[i+1])/(hr[i+1]+hr[i]);
					}

					//methods->getGradH(i+1,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k = 0; k < 2; k++)
					//	gradH[k] = (gradHL[k] * hr[i+1] + gradHR[k] * hr[i])/(hr[i+1] + hr[i]);

					//methods->Frvisc(Uaver,Frv2,tau,gradH[0]);
					methods->NMgdMethods::Roe(Ur,Ul,Fr2,"r","+");

					if(isnan(Fr2))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fr2 is nan!");
						exit(1);
					}

				}

				if(j == 0)
				{
					methods->get_grad1(i,j-1,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j-1,U_pack_arr,psiL,bz,br,hr,hz,"-");
					methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");

					//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j-1,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]+tauR[l][m])*.5;

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U1[k]+psiL[k]*gradUl[1][k]*hz[j]*.5;
						Ur[k]=U[i][j][k]-psiR[k]*gradUr[1][k]*hz[j]*.5;

						//Uaver[k]=(U1[k]+U[i][j][k])*.5;
					}

					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j-1,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]+gradHR[k])*.5;

					//methods->Fzvisc(Uaver,Fzv1,tau,gradH[1]);
					methods->NMgdMethods::Roe(Ur,Ul,Fz1,"z","-");

					if(isnan(Fz1))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz1 is nan!");
						exit(1);
					}

					methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j+1,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
					methods->limiter(i,j+1,U_pack_arr,psiR,bz,br,hr,hz,"+");

					//methods->stress(i,j+1,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hz[j+1]+tauR[l][m]*hz[j])/(hz[j+1]+hz[j]);

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j][k]+psiL[k]*gradUl[1][k]*hz[j]*.5;
						Ur[k]=U[i][j+1][k]-psiR[k]*gradUr[1][k]*hz[j+1]*.5;

						//Uaver[k]=(U[i][j+1][k]*hz[j]+U[i][j][k]*hz[j+1])/(hz[j+1]+hz[j]);
					}

					//methods->getGradH(i,j+1,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]*hz[j+1]+gradHR[k]*hz[j])/(hz[j+1]+hz[j]);

					//methods->Fzvisc(Uaver,Fzv2,tau,gradH[1]);
					methods->NMgdMethods::Roe(Ur,Ul,Fz2,"z","+");

					if(isnan(Fz2))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz2 is nan!");
						exit(1);
					}
				}
				else
				{
					methods->get_grad1(i,j-1,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j-1,U_pack_arr,psiL,bz,br,hr,hz,"-");
					methods->limiter(i,j,U_pack_arr,psiR,bz,br,hr,hz,"-");

					//methods->stress(i,j,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j-1,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]*hz[j]+tauR[l][m]*hz[j-1])/(hz[j-1]+hz[j]);

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j-1][k]+psiL[k]*gradUl[1][k]*hz[j-1]*.5;
						Ur[k]=U[i][j][k]-psiR[k]*gradUr[1][k]*hz[j]*.5;

						//Uaver[k]=(U[i][j-1][k]*hz[j]+U[i][j][k]*hz[j-1])/(hz[j-1]+hz[j]);
					}
						
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j-1,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]*hz[j]+gradHR[k]*hz[j-1])/(hz[j-1]+hz[j]);

					//methods->Fzvisc(Uaver,Fzv1,tau,gradH[1]);
					methods->NMgdMethods::Roe(Ur,Ul,Fz1,"z","-");

					if(isnan(Fz1))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz1 is nan!");
						exit(1);
					}

					methods->get_grad1(i,j,U_pack_arr,gradUl,bz,br,hr,hz);
					methods->get_grad1(i,j+1,U_pack_arr,gradUr,bz,br,hr,hz);

					methods->limiter(i,j,U_pack_arr,psiL,bz,br,hr,hz,"+");
					methods->limiter(i,j+1,U_pack_arr,psiR,bz,br,hr,hz,"+");

					//methods->stress(i,j+1,const_cast<const double***>(U),tauR,bz,br,hr,hz);
					//methods->stress(i,j,const_cast<const double***>(U),tauL,bz,br,hr,hz);

					//for(int l=0;l<3;l++)
					//	for(int m=0;m<3;m++)
					//		tau[l][m]=(tauL[l][m]+tauR[l][m])*.5;

					for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
					{
						Ul[k]=U[i][j][k]+psiL[k]*gradUl[1][k]*hz[j]*.5;
						Ur[k]=U1[k]-psiR[k]*gradUr[1][k]*hz[j]*.5;

						//Uaver[k]=(U1[k]+U[i][j][k])*.5;
					}
						
					//methods->getGradH(i,j+1,const_cast<const double***>(U),hr,hz,bz,br,gradHR);
					//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradHL);

					//for(int k=0;k<2;k++)
					//	gradH[k]=(gradHL[k]+gradHR[k])*.5;

					//methods->Fzvisc(Uaver,Fzv2,tau,gradH[1]);
					methods->NMgdMethods::Roe(Ur,Ul,Fz2,"z","+");

					if(isnan(Fz2))
					{
						NTracer::traceToTerminal(__FUNCTION__, __LINE__, "Fz2 is nan!");
						exit(1);
					}


				}
			delete []Ust;
			}
			}
			catch(int e)
			{
				cout << "i=" << i << "\n";
				cout << "j=" << j << "\n";
				exit(1);
			}

			dlr=hr[i];
			dlz=hz[j];
			methods->NMgdMethods::Sourse(U[i][j], S);
			//methods->getGradH(i,j,const_cast<const double***>(U),hr,hz,bz,br,gradH);
	  //      methods->stress(i,j,const_cast<const double***>(U),tau,bz,br,hr,hz);

			//methods->SourseVisc(U[i][j],Sv,tau,gradH[0]);
			
			for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
			{
				Re[i][j][k] = (Fr2[k] - Fr1[k])/dlr + (Fz2[k] - Fz1[k])/dlz + S[k]/rij;
				//Rv[i][j][k]=-(Frv1[k]*dlz-Frv2[k]*dlz)-(Fzv1[k]*dlr-Fzv2[k]*dlr)-Sv[k]*hr[i]*hz[j]/rij;
			}

			//double *traceData = new double[NInitial::getNcomp() + 2];
			//traceData[NInitial::getNcomp()] = i;
			//traceData[NInitial::getNcomp() + 1] = j;
			//for(int l = 0; l < NInitial::getNcomp(); l++)
			//	traceData[l] = Fr1[l];
			//NTracer::traceToFile(traceData, "double", NInitial::getNcomp() + 2);
			//for(int l = 0; l < NInitial::getNcomp(); l++)
			//	traceData[l] = Fr2[l];
			//NTracer::traceToFile(traceData, "double", NInitial::getNcomp() + 2);
			//for(int l = 0; l < NInitial::getNcomp(); l++)
			//	traceData[l] = Fz1[l];
			//NTracer::traceToFile(traceData, "double", NInitial::getNcomp() + 2);
			//for(int l = 0; l < NInitial::getNcomp(); l++)
			//	traceData[l] = Fz2[l];
			//NTracer::traceToFile(traceData, "double", NInitial::getNcomp() + 2);
			//delete []traceData;
		}

	delete []U1;
	delete []Ul;
	delete []Ur;
	delete []Fr1;
	delete []Fr2;
	delete []Fz1;
	delete []Fz2;
	//delete []S;
	for(size_t i = 0; i < 2; i++)
		delete []gradUr[i];
	delete []gradUr;
	for(size_t i = 0; i < 2; i++)
		delete []gradUl[i];
	delete []gradUl;

	delete []psiL;
	delete []Frv1;
	delete []Frv2;
	delete []Fzv1;
	delete []Fzv2;
	delete []psiR;
	//delete []Sv;
	delete U_pack;

	delete methods;
	delete bound;
}
void NCalcs::restrict(double ***U, double &t, const double *hr, const double *hz)
{
	double **tau;
	const double sigma1 = 0.28;

	tau = new double*[NInitial::get_xmax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		tau[i] = new double[NInitial::get_ymax()];

	NMgdMethods *methods = new NMgdMethods;

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
		{
			tau[i][j] = sigma1 * hr[i]/(abs(U[i][j][1]/U[i][j][0]) + methods->sound(U[i][j], "r"));
			tau[i][j] = min(tau[i][j], sigma1 * hz[j]/(abs(U[i][j][2]/U[i][j][0]) + methods->sound(U[i][j], "z")));
		}
	
	t = tau[0][0];

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
			t = min(t, tau[i][j]);

	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		delete []tau[i];
	delete []tau;

	delete methods;
}
void NCalcs::solve(double ***U, double &R2, double &R3, double &R4, double &R5, double &tau, const double *hr, const double *hz)
{
	double ***Re, Vol, ***Rv;

	Re = new double**[NInitial::get_xmax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		Re[i] = new double*[NInitial::get_ymax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			Re[i][j] = new double[NInitial::getNcomp()];

	Rv = new double**[NInitial::get_xmax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		Rv[i] = new double*[NInitial::get_ymax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			Rv[i][j] = new double[NInitial::getNcomp()];

	restrict(U, tau, hr, hz);

	resid(U, Re, hr, hz, Rv);
//Euler method

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
		{ 	
			for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				U[i][j][k] -= tau * Re[i][j][k];

			if(U[i][j][4] < 0)
			{
				cout << "Energy on cell i=" << i << ", j=" << j << ", is negative!";
				exit(1);
			}
		}



//H korrection
	//double *traceData = new double[NInitial::getNcomp() + 2];
	NMgdMethods *methods = new NMgdMethods;
	double ***korrH;

	korrH = methods->korrectH(const_cast<const double***>(U), hr, hz, tau);

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
		{
			U[i][j][4] += .5 * (pow(korrH[i][j][0], 2) + pow(korrH[i][j][1], 2) + pow(korrH[i][j][2], 2) - 
				pow(U[i][j][5], 2) - pow(U[i][j][6], 2) - pow(U[i][j][7], 2));
			for(int k = 5; k <= NInitial::getNcomp() - 1; k++)
				U[i][j][k] = korrH[i][j][k - 5];

			//traceData[NInitial::getNcomp()] = i;
			//traceData[NInitial::getNcomp() + 1] = j;
			//for(int l = 0; l < NInitial::getNcomp(); l++)
			//	traceData[l] = U[i][j][l];
			//NTracer::traceToFile(traceData, "double", NInitial::getNcomp() + 2);
			if(U[i][j][4] < 0)
			{
				cout << "Energy on cell i=" << i << ", j=" << j << ", is negative after H korrection!";
				exit(1);
			}

		}
	//delete []traceData;


//runge-kutta method
	//for(int i=0;i<=xmax-1;i++)
	//	for(int j=0;j<=ymax-1;j++)
	//	{ 	
	//		Vol=hz[j]*hr[i];
	//		for(int k=0;k<=Ncomp-1;k++)
	//			koeff[0][i][j][k] = -tau/Vol*(Re[i][j][k]+Rv[i][j][k]);
	//	}

		//double U1[xmax][ymax][Ncomp];
		//for(int i=0;i<=xmax-1;i++)
		//	for(int j=0;j<=ymax-1;j++)
		//		for(int k=0;k<=Ncomp-1;k++)
		//			U1[i][j][k] = U[i][j][k] + .5*koeff[0][i][j][k];

		//resid(U1, Re, hr, hz, Rv);

		//for(int i=0;i<=xmax-1;i++)
		//	for(int j=0;j<=ymax-1;j++)
		//	{ 	
		//		Vol=hz[j]*hr[i];
		//		for(int k=0;k<=Ncomp-1;k++)
		//			koeff[1][i][j][k] = -tau/Vol*(Re[i][j][k]+Rv[i][j][k]);
		//	}




		//for(int i=0;i<=xmax-1;i++)
		//	for(int j=0;j<=ymax-1;j++)
		//		for(int k=0;k<=Ncomp-1;k++)
		//			U1[i][j][k] = U[i][j][k] + .5*koeff[1][i][j][k];

		//resid(U1, Re, hr, hz, Rv);

		//for(int i=0;i<=xmax-1;i++)
		//	for(int j=0;j<=ymax-1;j++)
		//	{ 	
		//		Vol=hz[j]*hr[i];
		//		for(int k=0;k<=Ncomp-1;k++)
		//			koeff[2][i][j][k] = -tau/Vol*(Re[i][j][k]+Rv[i][j][k]);
		//	}




		//for(int i=0;i<=xmax-1;i++)
		//	for(int j=0;j<=ymax-1;j++)
		//		for(int k=0;k<=Ncomp-1;k++)
		//			U1[i][j][k] = U[i][j][k] + koeff[2][i][j][k];

		//resid(U1, Re, hr, hz, Rv);

		//for(int i=0;i<=xmax-1;i++)
		//	for(int j=0;j<=ymax-1;j++)
		//	{ 	
		//		Vol=hz[j]*hr[i];
		//		for(int k=0;k<=Ncomp-1;k++)
		//			koeff[3][i][j][k] = -tau/Vol*(Re[i][j][k]+Rv[i][j][k]);
		//	}

		//for(int i=0;i<=xmax-1;i++)
		//	for(int j=0;j<=ymax-1;j++)
		//	{ 	
		//		Vol=hz[j]*hr[i];
		//		for(int k=0;k<=Ncomp-1;k++)
		//			U[i][j][k] += 0.16666666666666666666666666666667*(koeff[0][i][j][k] + 2*koeff[1][i][j][k]
		//		+ 2*koeff[2][i][j][k] + koeff[3][i][j][k]);
		//	}


	R2=Re[0][0][0];

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
			for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				R2=max(R2,Re[i][j][k]);

	R3=Re[0][0][0];

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
			for(int k = 0; k <= NInitial::getNcomp() - 1; k++)
				R3=min(R3,Re[i][j][k]);

	//R4=Rv[0][0][0];

	//for(int i=0;i<=xmax-1;i++)
	//	for(int j=0;j<=ymax-1;j++)
	//		for(int k=0;k<=Ncomp-1;k++)
	//			R4=max(R4,Rv[i][j][k]);

	//R5=Rv[0][0][0];

	//for(int i=0;i<=xmax-1;i++)
	//	for(int j=0;j<=ymax-1;j++)
	//		for(int k=0;k<=Ncomp-1;k++)
	//			R5=min(R5,Rv[i][j][k]);

	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			delete []Re[i][j];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		delete []Re[i];
	delete []Re;

	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			delete []Rv[i][j];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		delete []Rv[i];
	delete []Rv;

	delete methods;

	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			delete []korrH[i][j];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		delete []korrH[i];
	delete []korrH;
}
