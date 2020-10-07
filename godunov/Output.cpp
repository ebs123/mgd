#include "includes\Output.h"

COutput::COutput(void)
{
}

COutput::~COutput(void)
{
}

void COutput::save2dPlot(const char* namefile, CMeshGenerator* mesh, double time, double** R, double** U, double** V, double** P, double** S, double** S_diff)
{
	FILE *file;

	file = fopen(namefile, "w");
	fprintf(file, "%s %lf %s\n", "TITLE = \"", time, "\"");
	if(S_diff != NULL)
	{
		fprintf(file, "%s\n", "variables= \"x\", \"y\", \"R\", \"U\", \"V\", \"P\", \"S\", \"S_diff\"");
	}
	else
	{
		fprintf(file, "%s\n", "variables= \"x\", \"y\", \"R\", \"U\", \"V\", \"P\", \"S\"");
	}
	fprintf(file, "%s %d %s %d %s\n", "ZONE T= \"\", I=", mesh->getNumCells()[1], ", J=", mesh->getNumCells()[0], ", F=POINT");

	double *x = mesh->getMeshComponent(0);
	double *y = mesh->getMeshComponent(1);

	for(int i = 0; i < mesh->getNumCells()[0]; i++)
		for(int j = 0; j < mesh->getNumCells()[1]; j++)
		{
			if(S_diff != NULL)
			{
				fprintf(file, "%lf %lf %lf %lf %lf %lf %lf %lf\n", x[i], y[j], R[i][j], U[i][j], V[i][j], P[i][j], S[i][j], S_diff[i][j]);
			}
			else
			{
				fprintf(file, "%lf %lf %lf %lf %lf %lf %lf\n", x[i], y[j], R[i][j], U[i][j], V[i][j], P[i][j], S[i][j]);
			}
		}

	fclose(file);
}

void COutput::save1dPlotXAxis(const char* namefile, CMeshGenerator* mesh, double time, int y_layer_num, double** R, double** U, double** V, double** P)
{
	FILE *file;

	file = fopen(namefile, "w");
	fprintf(file, "%s %lf %s\n", "TITLE = \"", time, "\"");
	fprintf(file, "%s\n", "variables= \"x\", \"R\", \"U\", \"P\", \"E\", \"S\"");

	double *x = mesh->getMeshComponent(0);
//s=log(P[i][y_layer_num] / pow(R[i][y_layer_num], GAMMA))
	for(int i = 0; i < mesh->getNumCells()[0]; i++)
	{
		fprintf(file, "%lf %lf %lf %lf %lf %lf\n", x[i], R[i][y_layer_num], U[i][y_layer_num], P[i][y_layer_num], .5 * R[i][y_layer_num] * (U[i][y_layer_num] * U[i][y_layer_num]) + P[i][y_layer_num]/(R[i][y_layer_num] * (GAMMA - 1)), 
			log(P[i][y_layer_num] / pow(R[i][y_layer_num], GAMMA)));
	}

	fclose(file);
}

void COutput::save1dPlotYAxis(const char* namefile, CMeshGenerator* mesh, double time, int x_layer_num, double** R, double** U, double** V, double** P)
{
	FILE *file;

	file = fopen(namefile, "w");
	fprintf(file, "%s %lf %s\n", "TITLE = \"", time, "\"");
	fprintf(file, "%s\n", "variables= \"x\", \"R\", \"U\", \"P\", \"E\", \"S\"");

	double *y = mesh->getMeshComponent(1);

	for(int j = 0; j < mesh->getNumCells()[1]; j++)
	{
		fprintf(file, "%lf %lf %lf %lf %lf %lf\n", y[j], R[x_layer_num][j], U[x_layer_num][j], P[x_layer_num][j], .5 * R[x_layer_num][j] * (U[x_layer_num][j] * U[x_layer_num][j]) + P[x_layer_num][j]/(R[x_layer_num][j] * (GAMMA - 1)), 
			log(P[x_layer_num][j] / pow(R[x_layer_num][j], GAMMA)));
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

void COutput::saveKineticEnergySpectrum(const char* namefile, double E, int k)
{
	FILE *file;
	file = fopen(namefile, "a+");

	fprintf(file, "%lf %d\n", E, k);

	fclose(file);
}

void COutput::saveData(const char* namefile, double data1, double data2)
{
	FILE *file;
	file = fopen(namefile, "a+");

	fprintf(file, "%lf %lf\n", data1, data2);

	fclose(file);
}