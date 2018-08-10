#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\NOther.h"
#include <sstream>
#else
#include "includes/NOther.h"
#endif

using namespace std;

int main(int argc, char **argv)
{
	NMethods *methods = new NMethods;
	NInitial *init = new NInitial;
	NOther *other = new NOther;
	NCalcs *numCalcs = new NCalcs;

	double ***U, *hz, *hr;

#if defined(WIN64) | defined(_WIN64) 
	init->loadInitialData("C:\\Documents and Settings\\config\\config.xml");
#elif defined(WIN32) | defined(_WIN32)
	init->loadInitialData("C:\\Documents and Settings\\config\\config.xml");
#endif

	U = new double**[NInitial::get_xmax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		U[i] = new double*[NInitial::get_ymax()];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			U[i][j] = new double[NInitial::getNcomp()];

	hz = new double[NInitial::get_ymax() + 3];
	hr = new double[NInitial::get_xmax() + 3];

	if(NInitial::getIndex() == 0)
	{
		methods->LoadGrid(hr, hz);
		init->initial(U, hr, hz);
	}
	else
	{
		methods->LoadGrid(hr, hz);
		other->load(U, "outF.dat");
	}

	/*vector<int> sizes(2);
	sizes[0] = NInitial::get_xmax();
	sizes[1] = NInitial::get_ymax();
	int dummy_num = 4;
	NArrPacker<double> *U_pack = new NArrPacker<double>(dummy_num, sizes, NInitial::getPeriodicName(), NInitial::getSlipName(), (void*)U);
	FILE *file;
	file = fopen("debug.dat", "w");
	double*** arr = (double***)U_pack->getPackArr();
	for(size_t i = 0; i < NInitial::get_xmax() + 2 * dummy_num; i++)
		for(size_t j = 0; j < NInitial::get_ymax() + 2 * dummy_num; j++)
		{
			for(size_t k = 0; k < NInitial::getNcomp(); k++)
				fprintf(file, "%f ", arr[i][j][k]);
			fprintf(file, "%d %d\n", i, j);
		}
	exit(1);*/
	NTracer::openFile("flowTrace.dat");

	//other->save(const_cast<const double***>(U), "out.dat", hr, hz);
    other->save2dr(const_cast<const double***>(U), "out2d.dat", 200, hr, hz);
	cout << "file is saved" << '\n';

	double R, R3, R4, R5, tau, t = 0;
	double *traceData = new double[10];

	for(int i = 1; i <= NInitial::getNstep(); i++)
	{
		numCalcs->solve(U, R, R3, R4, R5, tau, hr, hz);
		t = tau + t;
		cout << "Nstep = " << i << ", resmax = " << R << ", resmin = " << R3 << ", resmaxv = " << R4 << 
			", resminv = " << R5 << ", t = " << t << '\n';

		if((i % NInitial::getNsave()) == 0)
		{
			stringstream str;
			str << i << ".dat";
			other->save2dr(const_cast<const double***>(U), str.str().c_str(), 200, hr, hz);
			//other->save(const_cast<const double***>(U), str.str().c_str(), hr, hz);
			cout<<"file is saved"<<'\n';
		}
	}
	delete []traceData;

	NTracer::closeFile();
	delete methods;
	delete init;
	delete other;
	delete numCalcs;
	delete []hr;
	delete []hz;
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		for(size_t j = 0; j < NInitial::get_ymax(); j++)
			delete []U[i][j];
	for(size_t i = 0; i < NInitial::get_xmax(); i++)
		delete []U[i];
	delete []U;
}