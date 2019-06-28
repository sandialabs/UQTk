#ifndef TMCMCHEADERSEEN
#define TMCMCHEADERSEEN
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <functional>
#include <string>

#include "dsfmt_add.h"

#ifdef USE_HDF5
#include "H5Cpp.h"
#endif


typedef std::vector<double>  RealVector;
typedef std::vector<int>      IntVector;
typedef std::vector<bool>    BoolVector;
typedef std::vector<char>    CharVector;

double tmcmc(RealVector &rngs, double gm, int nspl, RealVector dtvec, int iseed,
             const int nProcs);

double tmcmc(RealVector &rngs, double gm, int nspl, RealVector dtvec, int iseed,
           const int nProcs, int ndim, bool spc, double cv, double a,
           double b, bool mala, double tauScale, double betaThres, int MFactor, bool basis, int CATSteps) ;

#endif
