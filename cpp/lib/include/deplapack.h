/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#ifndef LAPACK_H
#define LAPACK_H

#include "ftndefs.h"

// Computes eigenvalues and eigenvectors of a real symmetric matrix A
// See dsyevx.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(dsyevx)(char*, char*, char*, int*,
				      double*, int*, double*, double*,
				      int*, int*, double*, int*,
				      double*, double*, int*, double*, int*,
				      int*, int*, int*);

// Computes Cholesky decomposition of a real symmetric matrix A
// See dpotrf.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(dpotrf)(char*, int*, double*, int*, int*);

// Computes the inverse of a real symmetric matrix A
// See dpotri.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(dpotri)(char*, int*, double*, int*, int*);

// Computes the inverse of a complex matrix
// see zgetri.f for details on the arguments.
extern FTN_FUNC void FTN_NAME(zgetri)(int*, double*, int*, int*, double*, int*, int*);


// Computes the least squares solution of an over-determined Ax=b
// see dgels.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgels)(char*, int*, int*, int*, double*, int*, double*, int*, double*, int*, int*);

// Computes the LU factorization of a general matrix
// see dgetrf.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgetrf)(int*, int*, double*, int*, int*, int*);

// Computes the inverse of a general matrix using its LU factorization
// see dgetri.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgetri)(int*, double*, int*, int*, double*, int*, int*);

// Computes the solution to a real system of linear equations AX = B
// see dgesv.f for details on the arguments
extern FTN_FUNC void FTN_NAME(dgesv)(int*, int*, double*, int*, int*, double*, int*, int*);

//
extern FTN_FUNC void FTN_NAME(dstev)( char *, int *, double *, double *, double *, int *, double *, int* ) ;

//
extern FTN_FUNC void FTN_NAME(dsteqr)( char *, int *, double *, double *, double *, int *, double *, int* ) ;

/*****************************************
New routines
*****************************************/

//DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
extern FTN_FUNC void FTN_NAME(dgeqrf)(int*, int*, double*, int*, double*, double*, int*, int*);

// DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
extern FTN_FUNC void FTN_NAME(dorgqr)(int*, int*, int*, double*, int*, double*, double*, int*, int*);

// DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
//    $                   WORK, LWORK, INFO )
extern FTN_FUNC void FTN_NAME(dgesvd)(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);

// DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
extern FTN_FUNC void FTN_NAME(dposv)(char*, int*, int*, double*, int*, double*, int*, int*);


#endif  /* LAPACK_H */
