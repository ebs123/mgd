#pragma once
#include "NCalcs.h"

class NOther
{
public:
	NOther(void);
	~NOther(void);

	void save(const double ***U, const char *namefile, const double *hr, const double *hz);
	void convert_to_3D(const char *namefile_in, const char *namefile_out);
	void save2dr(const double ***U, const char *namefile, const int j, const double *hr, const double *hz);
	void save2dz(const double ***U, const char *namefile, const int i, const double *hr,
		const double *hz);
	void save_rot(const double &rot, FILE *file, const double &p_up, const double &p_dwn, const double &time);
	void load(double ***U, char *namefile);
	double rot_phi(const double ***U, char *parameter, const double *hr, const double *hz);
	double circulation(const double ***U, const double *hr, const double *hz);
	//double divH(const double ***U, const double *hr, const double *hz, const int &i, const int &j);

};

