#include "includes\Output.h"

COutput::COutput(void)
{
}

COutput::~COutput(void)
{
}

void COutput::save2dPlot(const char* namefile, CMeshGenerator* mesh, double time, double** R, double** U, double** V, double** P, double** S)
{
	FILE *file;

	file = fopen(namefile, "w");
	fprintf(file, "%s %lf %s\n", "TITLE = \"", time, "\"");
	fprintf(file, "%s\n", "variables= \"x\", \"y\", \"R\", \"U\", \"V\", \"P\", \"S\"");
	fprintf(file, "%s %d %s %d %s\n", "ZONE T= \"\", I=", mesh->getNumCells()[0], ", J=", mesh->getNumCells()[1], ", F=POINT");

	double *x = mesh->getMeshComponent(0);
	double *y = mesh->getMeshComponent(1);

	for(int i = 0; i < mesh->getNumCells()[0]; i++)
		for(int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			fprintf(file, "%lf %lf %lf %lf %lf %lf %lf\n", x[i], y[j], R[i][j], U[i][j], V[i][j], P[i][j], S[i][j]);
		}

	fclose(file);
}

void COutput::save1dPlotXAxis(const char* namefile, CMeshGenerator* mesh, double time, int y_layer_num, double** R, double** U, double** V, double** P, double** S)
{
	FILE *file;

	file = fopen(namefile, "w");
	fprintf(file, "%s %lf %s\n", "TITLE = \"", time, "\"");
	fprintf(file, "%s\n", "variables= \"x\", \"R\", \"U\", \"V\", \"P\", \"S\"");

	double *x = mesh->getMeshComponent(0);

	for(int i = 0; i < mesh->getNumCells()[0]; i++)
	{
		fprintf(file, "%lf %lf %lf %lf %lf %lf\n", x[i], R[i][y_layer_num], U[i][y_layer_num], V[i][y_layer_num], P[i][y_layer_num], S[i][y_layer_num]);
	}

	fclose(file);
}

void COutput::save1dPlotYAxis(const char* namefile, CMeshGenerator* mesh, double time, int x_layer_num, double** R, double** U, double** V, double** P, double** S)
{
	FILE *file;

	file = fopen(namefile, "w");
	fprintf(file, "%s %lf %s\n", "TITLE = \"", time, "\"");
	fprintf(file, "%s\n", "variables= \"y\", \"R\", \"U\", \"V\", \"P\", \"S\"");

	double *y = mesh->getMeshComponent(1);

	for(int j = 0; j < mesh->getNumCells()[1]; j++)
	{
		fprintf(file, "%lf %lf %lf %lf %lf %lf\n", y[j], R[x_layer_num][j], U[x_layer_num][j], V[x_layer_num][j], P[x_layer_num][j], S[x_layer_num][j]);
	}

	fclose(file);
}

void COutput::saveKineticEnergyAndEnstrophy(const char* namefile, CMeshGenerator* mesh, double** R, double** U, double** V, double time)
{
	FILE *file;
	double energy, enstrophy;

	energy = 0;
	enstrophy = 0;
	for(int i = 0; i < mesh->getNumCells()[0]; i++)
	{
		for(int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			energy += .5 * R[i][j] * (U[i][j] * U[i][j] + V[i][j] * V[i][j]) * mesh->getMeshStep(0) * 
				mesh->getMeshStep(1);
			if(i > 0 && i < mesh->getNumCells()[0] - 1 && j > 0 && j < mesh->getNumCells()[1] - 1)
			{
				enstrophy += mesh->getMeshStep(0) * mesh->getMeshStep(1) * pow((V[i + 1][j] - V[i - 1][j]) * .5 / mesh->getMeshStep(0) - 
					(U[i][j + 1] - U[i][j - 1]) * .5 / mesh->getMeshStep(1), 2);
			}
		}
	}

	file = fopen(namefile, "a+");
	fprintf(file, "%lf %lf %lf\n", energy, enstrophy, time);
	fclose(file);
}

void COutput::saveKineticEnergySpectrum(const char* namefile, double E_spectrum, double k_x, double k_y)
{
	FILE *file;
	file = fopen(namefile, "a+");

	fprintf(file, "%lf %lf\n", E_spectrum, sqrt(k_x * k_x + k_y * k_y));

	fclose(file);
}