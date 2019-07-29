%module(directors="1") bcs
//=====================================================================================
//                     The UQ Toolkit (UQTk) version @UQTKVERSION@
//                     Copyright (@UQTKYEAR@) Sandia Corporation
//                     http://www.sandia.gov/UQToolkit/
//
//     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
//     with Sandia Corporation, the U.S. Government retains certain rights in this software.
//
//     This file is part of The UQ Toolkit (UQTk)
//
//     UQTk is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Lesser General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     UQTk is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Lesser General Public License for more details.
//
//     You should have received a copy of the GNU Lesser General Public License
//     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
//
//     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
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




