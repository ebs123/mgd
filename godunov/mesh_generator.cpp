#include "includes\mesh_generator.h"

CMeshGenerator::mesh_generator(double *fp_domain_length, double *fp_mesh_size, int f_dimension) : mp_domain_length(fp_domain_length), 
	mp_mesh_size(fp_mesh_size), m_mesh_dimension(f_dimension)
{
	mp_x_i = new double*[m_mesh_dimension];
	for(size_t i = 0; i < m_mesh_dimension; i++)
	{
		mp_x_i[i] = new double[mp_mesh_size[i]];
	}
}

CMeshGenerator::~mesh_generator(void)
{
	delete []mp_domain_length;
	delete []mp_mesh_size;

	for(size_t i = 0; i < m_mesh_dimension; i++)
	{
		delete []mp_x_i[i];
	}
	delete []mp_x_i;
}

double* CMeshGenerator::getMeshComponent(int f_component_number)
{
	double dx = mp_domain_length[f_component_number]/(mp_mesh_size[f_component_number] - 1);

	for(size_t i = 0; i < mp_mesh_size[f_component_number]; i++)
	{
		mp_x_i[f_component_number][i] = (i + .5) * dx;
	}

	return mp_x_i[f_component_number];
}

int* CMeshGenerator::getNumCells()
{
	return mp_mesh_size;
}

double CMeshGenerator::getMeshStep(int f_component_number)
{
	return mp_domain_length[f_component_number]/(mp_mesh_size[f_component_number] - 1);
}