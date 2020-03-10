/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
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

double tmcmc (RealVector &spls, RealVector &lprior, RealVector &llik,
           double gm, int nspl, int iseed,
           int nProcs, int ndim, double cv, int MFactor,
           bool basis, int CATSteps, int write_flag) ;

#endif
