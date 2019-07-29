#ifndef SLATEC_H
#define SLATEC_H

#include "ftndefs.h"


// Some typedefs for the function callbacks in its arguments
typedef void (*f77_matvecprod)(int*, double*, double*, int*,
        int*, int*, double*, int*);
typedef void (*f77_precond)(int*, double*, double*, int*,
        int*, int*, double*, int*, double*, int*);

// Preconditioned GMRES iterative sparse Ax=b solver
// see dgmres.f for details
extern FTN_FUNC void FTN_NAME(dgmres)(int*, double*, double*, int*, int*, int*,
        double*, int*, f77_matvecprod, f77_precond,
        int*, double*, int*, int*, double*, int*, int*, double*, double*,
        double*, int*, int*, int*, double*, int*);

#endif /* SLATEC_H */
