#include "includes\Output.h"

COutput::COutput(void)
{
}

COutput::~COutput(void)
{
}

void COutput::save2dPlot(const char* namefile, CMeshGenerator* mesh, double** R, double** U, double** V, double** P, double** S)
{
	FILE *file;

	file = fopen(namefile, "w");
	fprintf(file, "%s\n", "TITLE = \" \"");
	fprintf(file, "%s\n", "variables= \"x\", \"y\", \"R\", \"U\", \"V\", \"P\", \"S\"");
	fprintf(file, "%s %d %s %d %s\n", "ZONE T=\" \", I=", mesh->getNumCells()[0], ", J=", mesh->getNumCells()[1], ", F=POINT");

	double *x = mesh->getMeshComponent(0);
	double *y = mesh->getMeshComponent(1);

	for(int i = 0; i < mesh->getNumCells()[0]; i++)
		for(int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			fprintf(file, "%lf %lf %lf %lf %lf %lf %lf\n", x[i], y[j], R[i][j], U[i][j], V[i][j], P[i][j], S[i][j]);
		}

	fclose(file);
}