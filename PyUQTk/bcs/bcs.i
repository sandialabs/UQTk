%module(directors="1") bcs
//=====================================================================================
//
//                      The UQ Toolkit (UQTk) version 3.1.1
//                          Copyright (2021) NTESS
//                        https://www.sandia.gov/UQToolkit/
//                        https://github.com/sandialabs/UQTk
//
//     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
//     retains certain rights in this software.
//
//     This file is part of The UQ Toolkit (UQTk)
//
//     UQTk is open source software: you can redistribute it and/or modify
//     it under the terms of BSD 3-Clause License
//
//     UQTk is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     BSD 3 Clause License for more details.
//
//     You should have received a copy of the BSD 3 Clause License
//     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
//
//     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
//     Sandia National Laboratories, Livermore, CA, USA
//=====================================================================================

%{
#define SWIG_FILE_WITH_INIT
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <math.h>
#include "../../cpp/lib/bcs/bcs.h"

%}

/*************************************************************
// Standard SWIG Templates
*************************************************************/

// Include standard SWIG templates
// Numpy array templates and wrapping
%include "pyabc.i"
%include "../numpy/numpy.i"
%include "std_vector.i"
%include "std_string.i"
%include "cpointer.i"

%init %{
    import_array();
%}

%pointer_functions(double, doublep);

/*************************************************************
// Numpy SWIG Interface files
*************************************************************/

/*************************************************************
// Include header files
*************************************************************/

%include "../../cpp/lib/bcs/bcs.h"


%include "bcs_ext.py"
