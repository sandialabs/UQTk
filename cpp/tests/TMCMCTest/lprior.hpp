#ifndef LPRIOR_H_SEEN
#define LPRIOR_H_SEEN
#include <iostream>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cassert>
#ifdef USE_HDF5
#include "H5Cpp.h"
#endif

typedef std::vector<int> IntVector;
typedef std::vector<double> RealVector;
typedef std::vector<char> CharVector;

double normPDF(double x, double mean, double var);
double unifPDF(double x, double a, double b);
double priorLL(RealVector &a, CharVector &distr, RealVector &paramA, RealVector &paramB);

void parseLP(std::string setup_fname, std::string sample_fname);

#endif
