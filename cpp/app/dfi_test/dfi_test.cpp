#include <iostream>
#include "math.h"
#include "Array1D.h"
#include "Array2D.h"
#include "mcmc.h"
#include "quad.h"
#include "dfi.h"
#include "dsfmt_add.h"
// #include </usr/local/include/mpi.h>
// #include <mpi.h> // for nersc
// #include </usr/local/mpich-g/include/mpi.h> // for grover

using namespace std; 

/*************************************************
Setup the initial data vectors: beta, z and sigma 
*************************************************/
void setup_data(Array1D<double>& beta0, Array1D<double>& z0, Array1D<double>& sigma0, Array1D<double>& zs0){

	int nSigma = 1;
	sigma0.Resize(nSigma);
	sigma0(0) = .5; 

	int nBeta = 2; 
	beta0.Resize(nBeta);
	beta0(0) = 3.2; beta0(1) = 1.1; 

	int nX = 25;
	z0.Resize(2*nX);
	double dx = 5.0/nX; 
	double x,y; 
	dsfmt_gv_init_gen_rand(13); // set initial seed to add noise to y
	for (int i = 0; i < nX; i++){
		x = i*dx;
		z0(i) = x;
		y = 1.1*x + 3.2;
		z0(i+nX) = y + sigma0(0)*dsfmt_gv_genrand_nrv(); 
	}
	// z0.DumpBinary("data/z0.dat");
	// beta0.DumpBinary("data/beta0.dat");
	// sigma0.DumpBinary("data/sigma0.dat");

	// create starting vector which combines (z,s)
	zs0.Resize(2*nX + nSigma);
	for (int i = 0; i < zs0.Length()-1; i++){
		zs0(i) = z0(i);
	} 
	zs0(2*nX + nSigma -1) = sigma0(0);
	// zs0.DumpBinary("data/zs0.dat");

}

class DFISetup: public DFISetupBase{
public:
	DFISetup(){};
	~DFISetup(){};

	double f(Array1D<double>&, Array1D<double>&, Array1D<double>&);
	double S(Array1D<MCMC::chainstate> inner_samples);
	double fun(double x){return 2*x;}

};
/*************************************************
Define the likelihood function f(beta,z,sigma) 
and the outer likelihood as a function of the 
inner samples, S(inner_samples)
*************************************************/
double DFISetup::f(Array1D<double>& beta, Array1D<double>& z, Array1D<double>& sigma){
	// unravel z to get x and y data from z = (x,y)
	int N = z.Length()/2;
	Array1D<double> x(N);
	Array1D<double> y(N);
	for (int i = 0; i < N; i++){
		x(i) = z(i);
		y(i) = z(N+i);
	}

	// unravel beta = (b,a) into intercept and slope terms
	double b = beta(0);
	double a = beta(1); 

	// unravel sigma
	double s = sigma(0);

	// determine squared error:
	// sum_i( (y_i - a*x_i - b)^2 ) 
	double sqerr = 0; 
	for (int i = 0; i < N; i++){
		sqerr += pow(y(i) - (a*x(i) + b),2);
	}

	// get flag if any x value is less than 0 or greater than 5
	int xflag = 0; 
	double buffer = .1; 
	for (int i = 0; i < N; i++){
		if (x(i) < 0 - buffer || x(i) > 5 + buffer){
			xflag = 1; 
			break;
		}
	}

	// return log posterior of linear regression model
	// if x is outside boundary, return -1e12 ~ -inf
	if (xflag == 1){
		return -1e16; 
	}
	else{
		return -(.5)*(1/s/s)*sqerr + (-N-3)*s;
	} 
}
double DFISetup::S(Array1D<MCMC::chainstate> inner_samples){

	// this function takes the inner chain samples and outputs 
	// some function of the summary statistics
	int nSigma = 1; 
	int nBeta = 2; 
	int nSamples = inner_samples.Length();
	int nBurn = 10000; // must be less than default 30000

	Array1D<double> means(nBeta,0.0);
	Array1D<double> stds(nBeta,0.0);

	for (int d = 0; d < nBeta; d++){
		for (int i = nBurn; i < nSamples; i++){
			means[d] += inner_samples(i).state[d];
			stds[d] += pow(inner_samples(i).state[d],2);
		}
		means[d] *= 1.0/(nSamples - nBurn);
		stds[d] *= 1.0/(nSamples - nBurn - 1);
		stds[d] += - pow(means[d],2);
		stds[d] = sqrt(stds[d]);
	}

	// set summary statistic, i.e. means and stnd deviation of beta
	Array1D<double> means0(nBeta);
	means0(0) = 3.2; means0(1) = 1.1; // set given means
	Array1D<double> stds0(nBeta);
	stds0(0) = .31/2; stds0(1) = .11/2; // set given standard deviation

	// set delta params
	double delta1 = 10, delta2 = 10; 

	// compute mean and variance check
	double d1 = 0, d2 = 0; 
	for (int i = 0; i < nBeta; i++){
		d1 += delta1*pow((means0(i) - means(i))/means0(i),2); // mean check
		d2 += delta2*pow((stds0(i) - stds(i))/stds0(i),2);  // std dev check
	} 

	return -(d1 + d2); 
}

/*************************************************
Begin main outer chain code
*************************************************/
int main(int argc, char ** argv){

	/*************************************************
	MPI Setup
	*************************************************/
	int mynode = 1, totalnodes;
	char fn[100];

	// MPI_Init(&argc,&argv); 
	// MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); 
	// MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

	/*************************************************
	Initial setup of beta, z, and sigma
	for linear regression model
	*************************************************/
	cout << "\nSetting up initial beta, z, and sigma variables..." << endl;

	// define dimensionality of beta, x and sigma 
	int nBeta = 2, nX = 25, nSigma = 1; 
	// initialize beta, z and sigma vectors
	Array1D<double> beta0, z0, sigma0, zs0;
	// run user-defined setup which defines starting points
	setup_data(beta0, z0, sigma0, zs0);

	/*************************************************
	Initiate DFI setup class
	User must define the function f(beta,z,sigma)
	and S(inner samples)
	*************************************************/
	cout << "\nInitializing DFI setup class which defines the log posterior and the summary statistic function" << endl;
	// Setup dfi
	DFISetup dfi_setup;

	// test log posterior for initial data
	cout << "Testing the logposterior function, f, of the initial variables:" << endl;
	cout << "\nf(beta0,z0,sigma0) is ";
	cout << dfi_setup.f(beta0,z0,sigma0) << "\n" << endl;


	// ************************************************
	// Initiate DFI inner class  and set inner start
	// *************************************************
	cout << "Passing the DFI setup object to the inner class, which we will obtain inner samples for beta" << "\n" << endl;
	DFIInner dfi_inner(dfi_setup);

	// set inner chain params (nCalls and nBurn are 30k and 10k by default)
	cout << "Setting inner chain params: dimension, init beta (prop width is adaptive)\n" << endl;
	dfi_inner.ndim = beta0.Length();
	dfi_inner.beta0_ = beta0; 

	// ********************************************
	// Now run outer loop over z and sigma
	// *********************************************
	// initiate DFI object
	DFI dfi(dfi_inner); 

	// set sigma dim and z dimension
	dfi.sdim = 1; 
	dfi.zdim = 2*nX;

	// define outer chain dimensions
	int ndim_outer = dfi.sdim + dfi.zdim;

	// define outer chain params
	int nCalls = 1; 
	Array1D<double> gammas(ndim_outer,.1);

	cout << "\n" << "Starting outer chain ..." << endl;
	dfi.runChain(nCalls, gammas, zs0, 13*(mynode+1), mynode);

	// MPI_Finalize();

	return 0; 

}
