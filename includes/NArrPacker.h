#pragma once
#include "bound.h"
/*
*	 ласс-упаковщик сетки в расширенную сетку с dummy-€чейками, 
*	индексы dummy-€чеек имеют значение [-n; -1] и [Nmax; Nmax + n - 1],
*	где Nmax - размерность сетки в данном направлении, n - количество dummy-€чеек
*/
template <class Type>
class NArrPacker
{
private:
	char *z_boundary_cond;
	char *r_boundary_cond;
	int dummy_num;//число dummy €чеек
	vector<int> dim_sizes;//размеры массива (сетки)
	void* packed_arr;//собственно, упакованный массив

public:
	NArrPacker(int num_dummy, vector<int> sizes, char *boundary_cond_z, char *boundary_cond_r, void* arr = NULL);
	virtual ~NArrPacker(void);

	void packer(void* arr);

	inline void* getPackArr();
};

template <class Type> void* NArrPacker<Type>::getPackArr()
{
	return packed_arr;
}

template <class Type> NArrPacker<Type>::NArrPacker(int num_dummy, vector<int> sizes, char *boundary_cond_z, char *boundary_cond_r, void* arr = NULL) : 
dummy_num(num_dummy), dim_sizes(sizes), z_boundary_cond(boundary_cond_z), r_boundary_cond(boundary_cond_r)
{
	if (sizes.size() == 1)
	{
		Type **tmp_arr;

		tmp_arr = new Type*[sizes[0] +  2 * dummy_num];
		for(size_t i = 0; i < sizes[0] + 2 * dummy_num; i++)
			tmp_arr[i] = new Type[NInitial::getNcomp()];

		packed_arr = (void*)tmp_arr;
	}
	else if(sizes.size() == 2)
	{
		Type ***tmp_arr;

		tmp_arr = new Type**[sizes[0] + 2 * dummy_num];
		for(size_t i = 0; i < sizes[0] + 2 * dummy_num; i++)
			tmp_arr[i] = new Type*[sizes[1] + 2 * dummy_num];
		for(size_t i = 0; i < sizes[0] + 2 * dummy_num; i++)
			for(size_t j = 0; j < sizes[1] + 2 * dummy_num; j++)
				tmp_arr[i][j] = new Type[NInitial::getNcomp()];

		packed_arr = (void*)tmp_arr;
	}
	else if(sizes.size() == 3)
	{
		Type ****tmp_arr;

		packed_arr = new Type***[sizes[0] + 2 * dummy_num];
		for(size_t i = 0; i < sizes[0] + 2 * dummy_num; i++)
			tmp_arr[i] = new Type**[sizes[1] + 2 * dummy_num];
		for(size_t i = 0; i < sizes[0] + 2 * dummy_num; i++)
			for(size_t j = 0; j < sizes[1] + 2 * dummy_num; j++)
				tmp_arr[i][j] = new Type*[sizes[2] + 2 * dummy_num];
		for(size_t i = 0; i < sizes[0] + 2 * dummy_num; i++)
			for(size_t j = 0; j < sizes[1] + 2 * dummy_num; j++)
				for(size_t k = 0; k < sizes[2] + 2 * dummy_num; k++)
					tmp_arr[i][j][k] = new Type[NInitial::getNcomp()];

		packed_arr = (void*)tmp_arr;
	}

	if(arr != NULL)
		packer(arr);

}

template <class Type> NArrPacker<Type>::~NArrPacker(void)
{
	if (dim_sizes.size() == 1)
	{
		Type **tmp_arr;
		tmp_arr = (Type**)packed_arr;

		for(size_t i = 0; i < dim_sizes[0] + 2 * dummy_num; i++)
			delete []tmp_arr[i];
		delete []tmp_arr;
	}
	else if(dim_sizes.size() == 2)
	{
		Type ***tmp_arr;
		tmp_arr = (Type***)packed_arr;

		for(size_t i = 0; i < dim_sizes[0] + 2 * dummy_num; i++)
			for(size_t j = 0; j < dim_sizes[1] + 2 * dummy_num; j++)
				delete []tmp_arr[i][j];
		for(size_t i = 0; i < dim_sizes[0] + 2 * dummy_num; i++)
			delete []tmp_arr[i];
		delete []tmp_arr;

	}
	else if(dim_sizes.size() == 2)
	{
		Type ****tmp_arr;
		tmp_arr = (Type****)packed_arr;

		for(size_t i = 0; i < dim_sizes[0] + 2 * dummy_num; i++)
			for(size_t j = 0; j < dim_sizes[1] + 2 * dummy_num; j++)
				for(size_t k = 0; k < dim_sizes[2] + 2 * dummy_num; k++)
					delete []tmp_arr[i][j][k];
		for(size_t i = 0; i < dim_sizes[0] + 2 * dummy_num; i++)
			for(size_t j = 0; j < dim_sizes[1] + 2 * dummy_num; j++)
				delete []tmp_arr[i][j];
		for(size_t i = 0; i < dim_sizes[0] + 2 * dummy_num; i++)
			delete []tmp_arr[i];

		delete []tmp_arr;
	}
}

template <class Type> void NArrPacker<Type>::packer(void* arr)
{
	if (dim_sizes.size() == 1)
	{
		Type **tmp_arr, **typed_arr;
		tmp_arr = (Type**)packed_arr;

		typed_arr = (Type**)arr;

		for(size_t i = 0; i < dim_sizes[0]; i++)
			for(size_t comp = 0; comp < NInitial::getNcomp(); comp++)
				tmp_arr[i + dummy_num][comp] = typed_arr[i][comp];

		//packed_arr = (void*)tmp_arr;

	}
	else if(dim_sizes.size() == 2)
	{
		Type ***tmp_arr, ***typed_arr;
		tmp_arr = (Type***)packed_arr;

		typed_arr = (Type***)arr;

		for(size_t i = 0; i < dim_sizes[0]; i++)
			for(size_t j = 0; j < dim_sizes[1]; j++)
				for(size_t comp = 0; comp < NInitial::getNcomp(); comp++)
					tmp_arr[i + dummy_num][j + dummy_num][comp] = typed_arr[i][j][comp];


		//заполнение значени€ми dummy €чеек
		NMgdBoundaryCond *bound = new NMgdBoundaryCond();
		if(!strcmp(z_boundary_cond, NInitial::getSlipName()))
		{
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t dummy = 1; dummy < dummy_num + 1; dummy++)
				{
					bound->slip(typed_arr[i][dummy_num - dummy], tmp_arr[i + dummy_num][dummy - 1], "z");
					bound->slip(typed_arr[i][NInitial::get_ymax() - 1 + dummy - dummy_num], 
						tmp_arr[i + dummy_num][NInitial::get_ymax() - dummy], "z");
				}

			for(size_t i = 1; i < dummy_num + 1; i++)
				for(size_t dummy = 1; dummy < dummy_num + 1; dummy++)
				{
					bound->slip(typed_arr[i - 1][dummy_num - dummy], tmp_arr[dummy_num - i][dummy - 1], "z");
					bound->slip(typed_arr[i - 1][NInitial::get_ymax() - 1 + dummy - dummy_num], 
						tmp_arr[dummy_num - i][NInitial::get_ymax() - dummy], "z");

					bound->slip(typed_arr[NInitial::get_xmax() - i][dummy_num - dummy], 
						tmp_arr[NInitial::get_xmax() + dummy_num + i - 1][dummy - 1], "z");
					bound->slip(typed_arr[NInitial::get_xmax() - i][NInitial::get_ymax() - 1 + dummy - dummy_num], 
						tmp_arr[NInitial::get_xmax() + dummy_num + i - 1][NInitial::get_ymax() - dummy], "z");
				}
		}
		else if(!strcmp(z_boundary_cond, NInitial::getPeriodicName()))
		{
			for(size_t i = 0; i < NInitial::get_xmax(); i++)
				for(size_t dummy = 1; dummy < dummy_num + 1; dummy++)
				{
					bound->periodic(i, dummy_num - dummy, typed_arr, tmp_arr[i + dummy_num][dummy - 1], "z");
					bound->periodic(i, NInitial::get_ymax() - 1 - dummy_num + dummy, typed_arr, 
						tmp_arr[i + dummy_num][NInitial::get_ymax() - dummy + 2 * dummy_num], "z");
				}

			//for(size_t i = 1; i < dummy_num + 1; i++)
			//	for(size_t dummy = 1; dummy < dummy_num + 1; dummy++)
			//	{
			//		bound->periodic(i - 1, dummy_num - dummy, typed_arr, tmp_arr[dummy_num - i][dummy - 1], "z");
			//		bound->periodic(i - 1, NInitial::get_ymax() - 1 - dummy_num + dummy, typed_arr, tmp_arr[dummy_num - i][NInitial::get_ymax() - dummy + 2 * dummy_num], "z");
			//	}
		}

		if(!strcmp(r_boundary_cond, NInitial::getSlipName()))
		{
			for(size_t j = 0; j < NInitial::get_ymax(); j++)
				for(size_t dummy = 1; dummy < dummy_num + 1; dummy++)
				{
					bound->slipinner(typed_arr[dummy_num - dummy][j], tmp_arr[dummy - 1][j + dummy_num]);
					bound->slip(typed_arr[NInitial::get_xmax() - 1 + dummy - dummy_num][j], 
						tmp_arr[NInitial::get_xmax() - dummy + 2 * dummy_num][j + dummy_num], "r");
				}

			for(size_t j = 1; j < dummy_num + 1; j++)
				for(size_t dummy = 1; dummy < dummy_num + 1; dummy++)
				{
					bound->slipinner(typed_arr[dummy_num - dummy][j - 1], tmp_arr[dummy - 1][dummy_num - j]);
					bound->slip(typed_arr[NInitial::get_xmax() - 1 + dummy - dummy_num][j - 1], 
						tmp_arr[NInitial::get_xmax() - dummy + 2 * dummy_num][dummy_num - j], "r");

					bound->slipinner(typed_arr[dummy_num - dummy][NInitial::get_ymax() - j], 
						tmp_arr[dummy - 1][NInitial::get_ymax() + dummy_num + j - 1]);
					bound->slip(typed_arr[NInitial::get_xmax() - 1 + dummy - dummy_num][NInitial::get_ymax() - j], 
						tmp_arr[NInitial::get_xmax() - dummy + 2 * dummy_num][NInitial::get_ymax() + dummy_num + j - 1], "r");
				}
		}
		else
		{
			for(size_t j = 0; j < NInitial::get_ymax(); j++)
				for(size_t dummy = 1; dummy < dummy_num + 1; dummy++)
				{
					bound->periodic(dummy - 1, j, typed_arr, tmp_arr[dummy - 1][j], "r");
					bound->periodic(NInitial::get_xmax() - dummy, j, typed_arr, 
						tmp_arr[dummy_num - 1 + NInitial::get_xmax() + dummy][j], "r");
				}
		}


		if(!strcmp(r_boundary_cond, NInitial::getPeriodicName()) & !strcmp(z_boundary_cond, NInitial::getPeriodicName()))
		{
			for(size_t i = 0; i < dummy_num; i++)
				for(size_t j = 0; j < dummy_num; j++)
				{
					bound->periodic(NInitial::get_xmax() + i, NInitial::get_ymax() + j, typed_arr, tmp_arr[NInitial::get_xmax() - i - 1][NInitial::get_ymax() - j - 1], "z");
					bound->periodic(NInitial::get_xmax() + i, dummy_num - j - 1, typed_arr, tmp_arr[NInitial::get_xmax() - i - 1][j], "z");
					bound->periodic(dummy_num - i - 1, dummy_num - j - 1, typed_arr, tmp_arr[i][j], "z");
					bound->periodic(dummy_num - i - 1, NInitial::get_ymax() + j, typed_arr, tmp_arr[i][NInitial::get_ymax() - j - 1], "z");
				}
		}


	}
	else if(dim_sizes.size() == 3)
	{
		Type ****tmp_arr, ****typed_arr;
		tmp_arr = (Type****)packed_arr;

		typed_arr = (Type****)arr;

		for(size_t i = 0; i < dim_sizes[0]; i++)
			for(size_t j = 0; j < dim_sizes[1]; j++)
				for(size_t k = 0; k < dim_sizes[2]; k++)
					for(size_t comp = 0; comp < NInitial::getNcomp(); comp++)
						tmp_arr[i + dummy_num][j + dummy_num][k + dummy_num][comp] = typed_arr[i][j][k][comp];
	}


}