#pragma once

class COutput
{
public:
	COutput(void);
public:
	~COutput(void);

	void save2dPlot(const char* namefile, CMeshGenerator* mesh, double** R, double** U, double** V, double** P, double** S);
};
