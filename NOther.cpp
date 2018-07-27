#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\NOther.h"
#else
#include "includes/NOther.h"
#endif

NOther::NOther(void)
{
}
NOther::~NOther(void)
{
}
void NOther::save(const double ***U, const char *namefile, const double *hr, const double *hz)
{
	FILE *file;
	NMethods *methods = new NMethods;

	double x, y, ux, uy, uz, p, ro;

	file = fopen(namefile, "w");
	fprintf(file, "%s\n", "TITLE = \" \"");
	fprintf(file, "%s\n", "variables= \"x\", \"y\", \"ux\", \"uy\", \"uz\", \"p\", \"ro\", \"Hx\", \"Hy\", \"Hz\"");
	fprintf(file, "%s %d %s %d %s\n", "ZONE T=\" \", I=", NInitial::get_ymax(), ", J=", NInitial::get_xmax(), ", F=POINT");
	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
		{
			x = methods->rij1(i, hr);
			y = methods->zij1(j, hz[0]);
			ux = U[i][j][1]/U[i][j][0];
			uy = U[i][j][2]/U[i][j][0];
			uz = U[i][j][3]/U[i][j][0];
			p = methods->press(U[i][j]);
			ro = U[i][j][0];
			fprintf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", x, y, ux, uy, uz, p, ro, U[i][j][5], U[i][j][6], U[i][j][7]);
		}

	fclose(file);

	delete methods;
}
void NOther::convert_to_3D(const char *namefile_in, const char *namefile_out)
{
	FILE *file_in, *file_out;
	const int n_K = 20;
	string str;
	ifstream ifs;
	int Ymax, Xmax;

	ifs.open(namefile_in);

	while(str[0] != 'I')
		ifs >> str;

	ifs >> str;
	Ymax = atoi(str.c_str());

	while(str[0] != 'J')
		ifs >> str;

	ifs >> str;
	Xmax = atoi(str.c_str());

	while(str.find_first_of(".") == string::npos)
		ifs >> str;

	file_out = fopen(namefile_out, "w");
	fprintf(file_out, "%s\n", "TITLE = \" \"");
	fprintf(file_out, "%s\n","variables= \"r\",\"y\",\"z\",\"ur\",\"uy\",\"uz\",\"p\",\"ro\"");
	fprintf(file_out, "%s %d %s %d %s %d %s\n", "ZONE T=\" \", I=", Ymax, ", J=", Xmax, ", K=", n_K,", F=POINT");
		
	double r, y, z, ur, uz, uphi, p, ro;

	while(ifs.good())
	{
		r = atof(str.c_str());
		ifs >> str;
		z = atof(str.c_str());
		ifs >> str;
		ur = atof(str.c_str());
		ifs >> str;
		uz = atof(str.c_str());
		ifs >> str;
		uphi = atof(str.c_str());
		ifs >> str;
		p = atof(str.c_str());
		ifs >> str;
		ro = atof(str.c_str());
		for(int k = 0; k < n_K; k++)
		{
			y = r*sin(k*1.5*3.14156/(n_K-1));
			double r1 = r*cos(k*1.5*3.14156/(n_K-1)), ux, uy;

			ux = ur*cos(k*1.5*3.14156/(n_K-1)) - uphi*sin(k*1.5*3.14156/(n_K-1));
			uy = ur*sin(k*1.5*3.14156/(n_K-1)) + uphi*cos(k*1.5*3.14156/(n_K-1));
				
			fprintf(file_out, "%lf %lf %lf %lf %lf %lf %lf %lf\n", r1, y, z, ux, uy, uz, p, ro);
		}
		ifs >> str;
	}



	fclose(file_out);
	ifs.close();
}
void NOther::save2dr(const double ***U, const char *namefile, const int j, const double *hr, const double *hz)
{
	FILE *file;
	NMethods *methods = new NMethods;

	double x, y, z, ux, uy, uz, uphi, p, ro;

	file = fopen(namefile, "w");
	fprintf(file, "%s\n", "TITLE = \" \"");
	fprintf(file, "%s\n", "variables= \"x\", \"ux\", \"uy\", \"uz\", \"p\", \"ro\", \"Hx\", \"Hy\", \"Hz\"");
	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
	{
		x = methods->rij1(i, hr);
		y = methods->zij1(j, hz[0]);
		ux = U[i][j][1]/U[i][j][0];
		uy = U[i][j][2]/U[i][j][0];
		uz = U[i][j][3]/U[i][j][0];
		p = methods->press(U[i][j]);
		ro  = U[i][j][0];
		fprintf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", x, ux, uy, uz, p, ro, U[i][j][5], U[i][j][6], U[i][j][7]);
	}

	fclose(file);
	delete methods;
}
void NOther::save2dz(const double ***U, const char *namefile, const int i, const double *hr, const double *hz)
{
	FILE *file;
	NMethods *methods = new NMethods;

	double r, z, ur, uz, uphi, p, ro;

	file=fopen(namefile, "w");
	fprintf(file, "%s\n", "TITLE = \" \"");
	fprintf(file, "%s\n", "variables= \"z\",\"ur\",\"uz\",\"uphi\",\"p\",\"ro\", \"Hr\", \"Hz\", \"Hphi\"");
	for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
	{
		r = methods->rij1(i, hr);
		z = methods->zij1(j, hz[0]);
		ur = U[i][j][1]/U[i][j][0];
		uz = U[i][j][2]/U[i][j][0];
		uphi = U[i][j][3]/U[i][j][0];
		p = methods->press(U[i][j]);
		ro = U[i][j][0];
		fprintf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", z, ur, uz, uphi, p, ro, U[i][j][5], U[i][j][6], U[i][j][7]);
	}
	fclose(file);
	delete methods;
}
void NOther::save_rot(const double &rot, FILE *file, const double &p_up, const double &p_dwn, const double &time)
{
	fprintf(file, "%f %f %f %f\n", rot, time, p_up, p_dwn);
}
void NOther::load(double ***U, char *namefile)
{
	FILE *file;
	char string1[30], string2[30], string3[30], string4[30], string5[30], string6[30], string7[30];
	double r, z, ur, uz, uphi, p, ro;

	if(!(file = fopen(namefile,"r")))
	{
		//cout<<namefile<<'\n';
		cout << "err in load file" << '\n';
		exit(1);
	}

	for(int i = 0; i <= NInitial::get_xmax() - 1; i++)
		for(int j = 0; j <= NInitial::get_ymax() - 1; j++)
		{
			fscanf(file,"%s %s %s %s %s %s %s\n", &string1, &string2, &string3, &string4, &string5,
				&string6, &string7);
			ur = atof(string3);
			uz = atof(string4);
			uphi = atof(string5);
			p = atof(string6);
			ro = atof(string7);
			U[i][j][0] = ro;
			U[i][j][1] = ro * ur;
			U[i][j][2] = ro * uz;
			U[i][j][3] = ro * uphi;
			U[i][j][4] = ro * (p/(NInitial::get_sigma() * ro)+.5 * (pow(ur, 2) + pow(uz, 2) + pow(uphi, 2)));
		}

	fclose(file);
}
double NOther::rot_phi(const double ***U, char *parameter, const double *hr, const double *hz)
{
	double ur_p, ur_m, uz_p, uz_m;
		
	if (parameter == "up")
	{
		ur_p = U[int(NInitial::get_xmax() * .5)][NInitial::get_ymax()][1]/U[int(NInitial::get_xmax() * .5)][NInitial::get_ymax()][0];
		ur_m = U[int(NInitial::get_xmax() * .5)][NInitial::get_ymax() - 1][1]/U[int(NInitial::get_xmax() * .5)][NInitial::get_ymax() - 1][0];
		uz_p = U[int(NInitial::get_xmax() * .5) + 1][NInitial::get_ymax()][1]/U[int(NInitial::get_xmax() * .5) + 1][NInitial::get_ymax()][0];
		ur_p = U[int(NInitial::get_xmax() * .5) - 1][NInitial::get_ymax()][1]/U[int(NInitial::get_xmax() * .5) - 1][NInitial::get_ymax()][0];
	}
	else if(parameter == "down")
	{
		ur_p = U[int(NInitial::get_xmax() * .5)][1][1]/U[int(NInitial::get_xmax() * .5)][1][0];
		ur_m = U[int(NInitial::get_xmax() * .5)][0][1]/U[int(NInitial::get_xmax() * .5)][0][0];
		uz_p = U[int(NInitial::get_xmax() * .5) + 1][0][1]/U[int(NInitial::get_xmax() * .5) + 1][0][0];
		ur_p = U[int(NInitial::get_xmax() * .5) - 1][0][1]/U[int(NInitial::get_xmax() * .5) - 1][0][0];
	}
	else
	{
		ur_p = U[int(NInitial::get_xmax() * .5)][int(NInitial::get_ymax() * .5) + 1][1]/U[int(NInitial::get_xmax() * .5)][int(NInitial::get_ymax() * .5) + 1][0];
		ur_m = U[int(NInitial::get_xmax() * .5)][int(NInitial::get_ymax() * .5) - 1][1]/U[int(NInitial::get_xmax() * .5)][int(NInitial::get_ymax() * .5) - 1][0];
		uz_p = U[int(NInitial::get_xmax() * .5) + 1][int(NInitial::get_ymax() * .5)][1]/U[int(NInitial::get_xmax() * .5) + 1][int(NInitial::get_ymax() * .5)][0];
		ur_p = U[int(NInitial::get_xmax() * .5) - 1][int(NInitial::get_ymax() * .5)][1]/U[int(NInitial::get_xmax() * .5) - 1][int(NInitial::get_ymax() * .5)][0];
	}

	return (ur_p - ur_m)/hz[NInitial::get_ymax()] - (uz_p - uz_m)/hr[NInitial::get_xmax()];
}
double NOther::circulation(const double ***U, const double *hr, const double *hz)
{
	double circ[4] = {0, 0, 0, 0};

	for(int j = 0; j < NInitial::get_ymax(); j++)
	{
		circ[0] = circ[0] + U[int(NInitial::get_xmax() * .5)][j][2]/U[int(NInitial::get_xmax() * .5)][j][0] * hz[10];//hz[10], because it's the homogenious mesh
		circ[2] = circ[2] + U[NInitial::get_xmax() - 1][j][2]/U[NInitial::get_xmax() - 1][j][0] * hz[10];//hz[10], because it's the homogenious mesh
	}

	for(int i = 0; i < NInitial::get_xmax(); i++)
	{
		circ[1] = circ[1] + U[i][0][1]/U[i][0][0] * hr[10];//hz[10], because it's the homogenious mesh
		circ[3] = circ[3] + U[i][NInitial::get_ymax() - 1][1]/U[i][NInitial::get_ymax() - 1][0] * hr[10];//hz[10], because it's the homogenious mesh
	}

	return circ[1] + circ[2] - circ[3] - circ[0];
}
/*double NOther::divH(const double ***U, const double *hr, const double *hz, const int &i, const int &j)
{
	//для равномерной сетки
	NMethods methods;
	NMgdBoundaryCond bound;

	if(i == )
	return .5/methods.rij1(i, hr) * (methods.rij1(i + 1, hr) * U[i + 1][j][5] - methods.rij1(i - 1, hr) * U[i - 1][j][5])/hr[i] + (U[i][j + 1][6] - U[i][j - 1][6])/hz[i];
}*/