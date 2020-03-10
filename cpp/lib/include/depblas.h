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
#ifndef BLAS_H
#define BLAS_H

#include "ftndefs.h"

// Computes the dot product for two vectors
// see ddot.f for details
extern FTN_FUNC double FTN_NAME(ddot)(int*, double*, int*, double*, int*);

// Copies the contents of a double* array into another double* array
// see dcopy.f for details
extern FTN_FUNC void FTN_NAME(dcopy)(int*, double*, int *, double*, int*);

// Scales a double* array by a constant
// see dscal.f for details
extern FTN_FUNC void FTN_NAME(dscal)(int*, double*, double*, int*);

// Computes constant times a vector plus a vector
// see daxpy.f for details
extern FTN_FUNC void FTN_NAME(daxpy)(int*, double*, double*, int*, double*, int*);

// Computes the Euclidean norm of a vector
// see dnrm2.f for details
extern FTN_FUNC double FTN_NAME(dnrm2)(int*, double*, int*);

// Symmetric matrix-Vector product
// see dsymv.f for details about the arguments
extern FTN_FUNC void FTN_NAME(dsymv)(char*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

// Performs linear matrix-vector operations
// see dgemv.f for details
extern FTN_FUNC void FTN_NAME(dgemv)(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

// Performs matrix-matrix operations
// see dgemm.f for details
extern FTN_FUNC void FTN_NAME(dgemm)(char*, char*, int*, int*, int*, double*, double*,  int*, double*, int*, double*, double*, int*);


#endif  /* BLAS_H */
