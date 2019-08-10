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


typedef std::vector<double>  RealVector;
typedef std::vector<int>     IntVector;
typedef std::vector<bool>    BoolVector;
typedef std::vector<char>    CharVector;

double tmcmc (std::vector<std::vector<double>> &rngs, double gm, int nspl, int iseed,
           int nProcs, int ndim, double cv, int MFactor,
           bool basis, int CATSteps) ;

#endif
