#ifndef LBFGS_H_Seen
#define LBFGS_H_Seen

#include "ftndefs.h"

/* This is a C library, not C++ */
#ifdef __cplusplus
extern "C"
{
#endif
int lbfgsDR(int n, int m, double *x, int *nbd, double *l, double *u,
            double (* func)(int, double *, void*),
            void   (* gradf)(int, double *, double *, void*), 
            void* userdata) ;
#ifdef __cplusplus
}
#endif

extern FTN_FUNC void FTN_NAME(setulb)(
           int *, int *, double *, double *, double *, int *, 
           double *, double *, double *, double *, double *, int *, 
           char *, int *, char *, unsigned int *, int *, double *);

#endif /* lbfgs */
