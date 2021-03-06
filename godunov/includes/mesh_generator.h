#pragma once

class CMeshGenerator
{
private:
	int m_mesh_dimension;
	double *mp_domain_length;
	int *mp_mesh_size;
	double **mp_x_i;

public:
	CMeshGenerator(double *fp_domain_length, int *fp_mesh_size, int f_dimension);

	~CMeshGenerator(void);

	double* getMeshComponent(int f_component_number);
	double getMeshStep(int f_component_number);
	int* getNumCells();
};
