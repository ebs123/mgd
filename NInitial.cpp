#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\NInitial.h"
#else
#include "includes/NInitial.h"
#endif

int NInitial::Nsave, NInitial::Nstep, NInitial::Index, NInitial::RZONE, NInitial::Ncomp, 
	NInitial::xmax, NInitial::ymax;
double NInitial::Cp, NInitial::kappa, NInitial::RFAKTOR, NInitial::roin, NInitial::uphiin, 
	NInitial::uzin, NInitial::urin, NInitial::pin, NInitial::sigma, NInitial::deltaz, NInitial::deltar, 
	NInitial::R1, NInitial::mu, NInitial::H[3];
multimap<string, string> NInitial::xmlTagOptionsMap;
map<string, double> NInitial::initialData;
vector<string> NInitial::tagsList;
const double NInitial::pi = 3.1415926535897932384626433832795;
const char *NInitial::slip_name = "slip", *NInitial::periodic_name = "per";

NInitial::NInitial(void)
{
	xmlTagOptionsMap.insert(pair<string, string>("<grid", "xmax"));
	xmlTagOptionsMap.insert(pair<string, string>("<grid", "ymax"));
	xmlTagOptionsMap.insert(pair<string, string>("<grid", "RFAKTOR"));
	xmlTagOptionsMap.insert(pair<string, string>("<grid", "R1"));
	xmlTagOptionsMap.insert(pair<string, string>("<grid", "deltar"));
	xmlTagOptionsMap.insert(pair<string, string>("<grid", "deltaz"));
	xmlTagOptionsMap.insert(pair<string, string>("<grid", "RZONE"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "roin"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "uphiin"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "uzin"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "urin"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "pin"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "Hr"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "Hz"));
	xmlTagOptionsMap.insert(pair<string, string>("<inidata", "Hphi"));
	xmlTagOptionsMap.insert(pair<string, string>("<flowSettings", "Cp"));
	xmlTagOptionsMap.insert(pair<string, string>("<flowSettings", "kappa"));
	xmlTagOptionsMap.insert(pair<string, string>("<flowSettings", "sigma"));
	xmlTagOptionsMap.insert(pair<string, string>("<flowSettings", "mu"));
	xmlTagOptionsMap.insert(pair<string, string>("<options", "Nsave"));
	xmlTagOptionsMap.insert(pair<string, string>("<options", "Nstep"));
	xmlTagOptionsMap.insert(pair<string, string>("<options", "Index"));
	xmlTagOptionsMap.insert(pair<string, string>("<options", "Ncomp"));

	tagsList.push_back("<grid");
	tagsList.push_back("<inidata");
	tagsList.push_back("<flowSettings");
	tagsList.push_back("<options");
}
NInitial::~NInitial(void)
{
}
double NInitial::domega(const double &r)//,const double R)
{
	const double domegamax = .1;
	double R2, alpha/*1,alpha2*/;

	R2 = R1 + deltar;
	//alpha1=domegamax/deltar;
	alpha = domegamax;

	if(r >= (R1 + R2) * .5)
		return -alpha*(r - R2);
			//return alpha*(r - .5*(R1 + R2));
	else
		return alpha*(r - R1);
			//return alpha*(.5*(R1 + R2) - r);
}
double NInitial::ad(const double &r)
{
	const double deltap = 3;
	double alpha, R2/*, R3, R4*/ /*alpha2*//*delta*/;

	R2 = R1 + deltar;
		//R3=R1+.4*deltar;
		//R4=R2-.4*deltar;
		//delta = (R2 - R1)/2;
		//alpha1=1.0/((R1 + R2)*.5 - R3);
		//alpha = 1/(R4 - R3);

	return deltap * .25/mu * ((pow(R2, 2) - pow(R1 ,2))*log(R2/r)/log(R2/R1) + pow(r, 2) - pow(R2, 2));

		//if((r > R1) && (r <= R3))
		//	return -1;
		//else if((r > R3) & (r <= R4))
		//	return -alpha*(R4 - r);
		//else
		//	return 0;

		//if((r > R1) && (r <= (R1 + R2)*.5))
		//	return 1.;
		//else
		//	return 0.;
}
double NInitial::adr(const double &r, const double &z)
{
	const double drmax = .01, pi = 3.1415926535897932384626433832795, n = 11;
	double R2, alpha, delta;

	R2 = R1 + deltar;
	delta = deltar/xmax * 10.;
	alpha = drmax/delta;

	if((r >= (R1 + R2) * .5) & (r < (R1 + R2)*.5 + delta))
		return -alpha * (r - ((R1 + R2)*.5 + delta))*sin(n*pi*z);
			//return alpha*(r - .5*(R1 + R2));
	else if((r < (R1 + R2)*.5) & (r > (R1 + R2)*.5 - delta))
		return alpha * (r - ((R1 + R2) * .5 - delta))*sin(n*pi*z);
	else
		return 0;
			//return alpha*(.5*(R1 + R2) - r);
}
double NInitial::U_phi(const double &z, const double &r, double &omega)
{
	const double delta = deltaz * .01, U_phi_max = 1., Omega = 1.;
	double alpha = U_phi_max/delta;

	if(z < delta)
	{
		omega = alpha * z * Omega;
		return alpha * z * Omega * r;
	}

	omega = U_phi_max * Omega;
	return U_phi_max * Omega * r;
}
void NInitial::initial(double ***U, const double *hr, const double *hz)
{
		//const double sigma=.4;
		//unsigned int k=0,l=90;
	NMethods *methods = new NMethods;

	if(methods == NULL)
		cout << "Can't create methods";
	
	double rij, uz, ur, uphi, zij, p;
	/****************OrszagЦTang MHD turbulence problem
	//const double gamma = 1.6666666666666666666666666666667;
	/********************************************
	/********** распад ћ√ƒ-разрыва ******************/
	const double roLeft = 1, roRight = .125, pLeft = 1, pRight = .1, HzLeft = 1, HzRight = -1, Hr = .75;
	/***********************************************/
	/****************Kelvin-Helmgoltz instability
	const double a = 1, u0 = 2, u0_tilda = .008, lambda = 15.707963267948966192313216916398;
	**********************************************/

	for(int i = 0; i <= xmax - 1; i++)
		for(int j = 0; j <= ymax - 1; j++)
		{
			//R2 = R1 + deltar;
			rij = methods->rij1(i, hr);
			zij = methods->zij1(j, hz[0]);
			/*****************OrszagЦTang MHD turbulence problem*/
			/*U[i][j][0] = pow(gamma, 2);
			U[i][j][1] = - U[i][j][0] * sin(zij);
			U[i][j][2] = U[i][j][0] * sin(rij);
			U[i][j][3] = 0.;
			U[i][j][5] = - sin(zij);
			U[i][j][6] = sin(2 * rij);
			U[i][j][7] = 0.;
			U[i][j][4] = gamma/sigma + (pow(U[i][j][1]/U[i][j][0], 2) + pow(U[i][j][2]/U[i][j][0], 2)) * 
				.5 * U[i][j][0] + .5 * (pow(U[i][j][5], 2) + pow(U[i][j][6], 2));*/
			
			/******************/
			/*****************Kelvin-Helmgoltz instability (нужно допилить граничные услови€)
			U[i][j][0] = 1.;
			U[i][j][1] = U[i][j][0] * .5 * u0 * tanh(zij/a);
			U[i][j][2] = 0.;
			U[i][j][3] = 0.;
			U[i][j][4] = p0/sigma + pow(.5 * u0 * tanh(zij/a), 2) * .5 * U[i][j][0] + .5 * pow(Hz0, 2);
			U[i][j][5] = 0.;
			U[i][j][6] = 0.;
			U[i][j][7] = Hz0;
			*****************/
			/********** распад ћ√ƒ-разрыва ******************/

			uz = 0.;
			ur = 0.;
			uphi = 0.;
			if(rij < 1)
			{
				p = pLeft;
				U[i][j][6] = HzLeft;
				U[i][j][0]=roLeft;
				U[i][j][1]=roLeft*ur;
				U[i][j][2]=roLeft*uz;
				U[i][j][3]=roLeft*uphi;
				U[i][j][4] = p/sigma + (pow(ur, 2) + pow(uz, 2) + pow(uphi, 2)) * .5 * roLeft + .5 * (pow(HzLeft, 2) + pow(H[0], 2) + pow(Hr, 2));
			}
			else
			{
				p=pRight;
				U[i][j][6] = HzRight;
				U[i][j][0]=roRight;
				U[i][j][1]=roRight*ur;
				U[i][j][2]=roRight*uz;
				U[i][j][3]=roRight*uphi;
				U[i][j][4]=p/sigma + (pow(ur, 2) + pow(uz, 2) + pow(uphi, 2)) * .5 * roRight + .5 * (pow(HzRight, 2) + pow(H[0], 2) + pow(Hr, 2));
			}

			U[i][j][7] = H[0];
			U[i][j][5] = Hr;
			/******************************/
		}

	delete methods;
}
int NInitial::getNsave()
{
	return Nsave;
}
int NInitial::getNstep()
{
	return Nstep;
}
int NInitial::getIndex()
{
	return Index;
}
int NInitial::getRZONE()
{
	return RZONE;
}
int NInitial::getNcomp()
{
	return Ncomp;
}
int NInitial::get_xmax()
{
	return xmax;
}
int NInitial::get_ymax()
{
	return ymax;
}
double NInitial::getCp()
{
	return Cp;
}
double NInitial::get_kappa()
{
	return kappa;
}
double NInitial::getRFAKTOR()
{
	return RFAKTOR;
}
double NInitial::get_roin()
{
	return roin;
}
double NInitial::get_uphiin()
{
	return uphiin;
}
double NInitial::get_uzin()
{
	return uzin;
}
double NInitial::get_urin()
{
	return urin;
}
double NInitial::get_pin()
{
	return pin;
}
double NInitial::get_sigma()
{
	return sigma;
}
double NInitial::get_deltaz()
{
	return deltaz;
}
double NInitial::get_deltar()
{
	return deltar;
}
double NInitial::getR1()
{
	return R1;
}
double NInitial::get_mu()
{
	return mu;
}
void NInitial::loadInitialData(const string &patch)
{
	fstream fileStream;
	char *pCharStr;

	pCharStr = new char[300];

	fileStream.open(patch.c_str(), fstream::in);

	if(!fileStream.good())
	{
		printf("Can't open config file: ");
		if(fileStream.eof())
			printf("eof\n");
		else if(fileStream.bad())
			printf("bad\n");
		else if(fileStream.fail())
			printf("fail\n");
		exit(0);
	}

	while(!fileStream.eof())
	{
		fileStream.getline(pCharStr, 300);
		string str(pCharStr);
		parseXMLTag(str);
	}

	for(map<string, double>::iterator it = initialData.begin(); it != initialData.end(); it++)
	{
		if(!strcmp((it->first).c_str(), "<grid xmax"))
		{
			xmax = it->second;
			cout << "xmax = " << xmax << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<grid ymax"))
		{
			ymax = it->second;
			cout << "ymax = " << ymax << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<grid RFAKTOR"))
		{
			RFAKTOR = it->second;
			cout << "RFAKTOR = " << RFAKTOR << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<grid R1"))
		{
			R1 = it->second;
			cout << "R1 = " << R1 << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<grid deltar"))
		{
			deltar = it->second;
			cout << "deltar = " << deltar << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<grid deltaz"))
		{
			deltaz = it->second;
			cout << "deltaz = " << deltaz << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<grid RZONE"))
		{
			RZONE = it->second;
			cout << "RZONE = " << RZONE << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata roin"))
		{
			roin = it->second;
			cout << "roin = " << roin << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata uphiin"))
		{
			uphiin = it->second;
			cout << "uphiin = " << uphiin << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata uzin"))
		{
			uzin = it->second;
			cout << "uzin = " << uzin << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata urin"))
		{
			urin = it->second;
			cout << "urin = " << urin << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata pin"))
		{
			pin = it->second;
			cout << "pin = " << pin << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata Hr"))
		{
			H[0] = it->second;
			cout << "Hr = " << H[0] << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata Hz"))
		{
			H[1] = it->second;
			cout << "Hz = " << H[1] << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<inidata Hphi"))
		{
			H[2] = it->second;
			cout << "Hphi = " << H[2] << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<flowSettings Cp"))
		{
			Cp = it->second;
			cout << "Cp = " << Cp << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<flowSettings kappa"))
		{
			kappa = it->second;
			cout << "kappa = " << kappa << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<flowSettings sigma"))
		{
			sigma = it->second;
			cout << "sigma = " << sigma << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<flowSettings mu"))
		{
			mu = it->second;
			cout << "mu = " << mu << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<options Nsave"))
		{
			Nsave = it->second;
			cout << "Nsave = " << Nsave << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<options Nstep"))
		{
			Nstep = it->second;
			cout << "Nstep = " << Nstep << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<options Index"))
		{
			Index = it->second;
			cout << "Index = " << Index << '\n';
		}
		else if(!strcmp((it->first).c_str(), "<options Ncomp"))
		{
			Ncomp = it->second;
			cout << "Ncomp = " << Ncomp << '\n';
		}

	}

	fileStream.close();
	delete[] pCharStr;
}
int NInitial::parseXMLTag(string &str)
{
	pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itMultimapPair;

	for(size_t l = 0; l < tagsList.size(); l++)
	{
		itMultimapPair = xmlTagOptionsMap.equal_range(tagsList[l]);
		for(multimap<string, string>::iterator itMultimap = itMultimapPair.first; itMultimap != itMultimapPair.second; itMultimap++)
		{
			if(str.find(itMultimap->first) != string::npos)
			{
				int pos;

				for(size_t i = 0; i < xmlTagOptionsMap.size(); i++)
				{
					if((pos = str.find(itMultimap->second)) == string::npos)
						printf("Error in cfg parsing: can't find %s %s tag option", (itMultimap->first).c_str(), (itMultimap->second).c_str());
					else
					{
						double value;
						pos = str.find('=', pos);
						char *charValue;

						if(str.find(' ') == pos + 1)
						{
							int charValueSize = str.find(' ', pos + 2) - pos - 2;
							charValue = new char[charValueSize];
							str.copy(charValue, charValueSize, pos + 2);
						}
						else
						{
							int charValueSize = str.find(' ', pos) - pos - 1;
							charValue = new char[charValueSize];
							str.copy(charValue, charValueSize, pos + 1);
						}
			#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
							value = atof(charValue);
			#else
							stringstream sstr;
							sstr << charValue;
							sstr >> value;
			#endif
							initialData.insert(pair<string, double>(itMultimap->first + ' ' + itMultimap->second, value));
							delete []charValue;
					}
				}
			}
		}
	}

	return 0;
}
double NInitial::get_pi()
{
	return pi;
}
double* NInitial::getH()
{
	return H;
}
char* NInitial::getPeriodicName()
{
	return const_cast<char*>(periodic_name);
}
char* NInitial::getSlipName()
{
	return const_cast<char*>(slip_name);
}