%module(directors="1") uqtkarray
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

%feature("autodoc", "3");
%rename(Assign) *::operator=;
%ignore *::operator[];

%{
#define SWIG_FILE_WITH_INIT
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <math.h>
#include "../../cpp/lib/array/Array1D.h"
#include "../../cpp/lib/array/Array2D.h"
#include "../../cpp/lib/array/arrayio.h"
#include "../../cpp/lib/array/arraytools.h"
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

// Basic typemap for an Arrays and its length.
// Must come before %include statement below

// For Array1D setnumpyarray4py function
%apply (long* IN_ARRAY1, int DIM1) {(long* inarray, int n)}
%apply (double* IN_ARRAY1, int DIM1) {(double* inarray, int n)}
// get numpy int and double array
%apply (long* INPLACE_ARRAY1, int DIM1) {(long* outarray, int n)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* outarray, int n)}

// For Array2D numpysetarray4py function
%apply (double* IN_FARRAY2, int DIM1, int DIM2) {(double* inarray, int n1, int n2)}
// get numpy array (must be FARRAY)
%apply (double* INPLACE_FARRAY2, int DIM1, int DIM2) {(double* outarray, int n1, int n2)}
// For Array2D numpysetarray4py function
%apply (long* IN_FARRAY2, int DIM1, int DIM2) {(long* inarray, int n1, int n2)}
// get numpy array (must be FARRAY)
%apply (long* INPLACE_FARRAY2, int DIM1, int DIM2) {(long* outarray, int n1, int n2)}


// For mcmc test to get log probabilities
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* l, int n)}

/*************************************************************
// Include header files
*************************************************************/

// // The above typemap is applied to header files below
%include "../../cpp/lib/array/Array1D.h"
%include "../../cpp/lib/array/Array2D.h"
%include "../../cpp/lib/array/arrayio.h"
%include "../../cpp/lib/array/arraytools.h"

// Typemaps for standard vector
// Needed to prevent to memory leak due to lack of destructor
// must use namespace std
namespace std{
    %template(dblVector) vector<double>;
    %template(intVector) vector<int>;
    %template(strVector) vector<string>;

}

%template(subMatrix_row_int) subMatrix_row<int>;
%template(subMatrix_row_dbl) subMatrix_row<double>;

%include "arrayext.i"

%pythoncode %{
import numpy as np
def uqtk2numpy(x):
	if x.type() == 'int':
		s = x.shape()
		imin = np.argmin(s)
		if len(s) == 1:
			n = s[0]
			y = np.zeros(n,dtype='int64')
			x.getnpintArray(y)
		if len(s) == 2 and np.amin(s) > 1:
			n = s[0]
			m = s[1]
			y = np.zeros((n,m),dtype='int64')
			x.getnpintArray(y)
		if len(s) == 2 and np.amin(s) == 1:
			y = np.array(x.flatten())
		return y.copy()
	else:
		s = x.shape()
		imin = np.argmin(s)
		if len(s) == 1:
			n = s[0]
			y = np.zeros(n)
			x.getnpdblArray(y)
		if len(s) == 2 and np.amin(s) > 1:
			n = s[0]
			m = s[1]
			y = np.zeros((n,m))
			x.getnpdblArray(y)
		if len(s) == 2 and np.amin(s) == 1:
			y = np.array(x.flatten())
		return y.copy()

def numpy2uqtk(y):
	s = np.shape(y)
	if len(s) == 1:
		n = s[0]
		x = dblArray1D(n)
	if len(s) == 2:
		n = s[0]
		m = s[1]
		x = dblArray2D(n,m)
	x.setnpdblArray(np.asfortranarray(y.copy()))
	return x
%}
