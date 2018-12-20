#pragma once
#include "mesh_generator.h"
#include <stdio.h>

class COutput
{
public:
	COutput(void);
public:
	~COutput(void);

	void save2dPlot(const char* namefile, CMeshGenerator* mesh, double** R, double** U, double** V, double** P, double** S);
	void save1dPlotXAxis(const char* namefile, CMeshGenerator* mesh, int y_layer_num, double** R, double** U, double** V, double** P, double** S);
	void save1dPlotYAxis(const char* namefile, CMeshGenerator* mesh, int x_layer_num, double** R, double** U, double** V, double** P, double** S);
};
