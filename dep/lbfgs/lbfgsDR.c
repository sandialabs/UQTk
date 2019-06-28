/*
    interface for lbfgs library.

    customized from driver2.f code downloaded with the library

     References:

        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
        memory algorithm for bound constrained optimization'',
        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
        Subroutines for Large Scale Bound Constrained Optimization''
        Tech. Report, NAM-11, EECS Department, Northwestern University,
        1994.

 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "lbfgs_routines.h"
void macheps() ;
void gradfnum(int n, double *x, double *xp, double (* func)(int, double *, void*), 
              double *g, void *userdata) ;

#define CSIZE     60
#define GRADSMALL 1.e-15
#define MAXFEVAL  100000

static double lbfgs_rt, lbfgs_at ;

int lbfgsDR(int n, int m, double *x, int *nbd, double *l, double *u,
            double (* func)(int, double *, void*),
            void   (* gradf)(int, double *, double *, void*), 
            void* userdata) 
{

  char task[CSIZE],csave[CSIZE] ;
  unsigned int lsave[4] ;

  int iprint, *iwa, isave[44], lenwa ;
  double f, factr, pgtol, *g, *xp, dsave[29], *wa ;

  int    i, lcont ;

  int ans = 0 ;

  lenwa = 2*m*n + 5*n + 12*m*m + 12*m ;

  iwa = (int    *)malloc(3*n*  sizeof(int)   ) ;
  wa  = (double *)malloc(lenwa*sizeof(double)) ;
  g   = (double *)malloc(n    *sizeof(double)) ;
  xp  = (double *)malloc(n    *sizeof(double)) ;
  
   
  /* We suppress the default output. */
  iprint = -1 ;

  /* 
     We suppress both code-supplied stopping tests because the
     user is providing his own stopping criteria. 
  */
  factr = 0.0 ;
  pgtol = 0.0 ;

  /* We start the iteration by initializing task. */
  memset(task,32,CSIZE) ; /* put spaces (ascii code 32 in task) */
  sprintf(task,"START") ;

  /* Compute machine epsilon */
  macheps() ;
  printf("Machine eps rt=%le, at=%le\n",lbfgs_rt,lbfgs_at);

  /* the beginning of the loop */
  do 
  {

    lcont = 0 ;
    setulb_(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task,
            &iprint, csave, lsave, isave, dsave) ;
    //printf("%s\n",task) ;
    if ( strncmp(task,"FG",2) == 0 )
    {
      
      /*
        the minimization routine has returned to request the
        function f and gradient g values at the current x.
      */

      /* Compute function value f for the sample problem. */
      f = func(n,x,userdata) ;

      /* Compute gradient g for the sample problem. */
      if ( gradf == NULL ){
        gradfnum(n,x,xp,func,g,userdata) ;
      }
      else{
        gradf(n, x, g, userdata) ;
      }
      /* go back to the minimization routine. */
      lcont = 1 ;
    }
    if ( strncmp(task,"NEW_X",5) == 0 )
    {
      /*
        the minimization routine has returned with a new iterate.
        At this point have the opportunity of stopping the iteration 
        or observing the values of certain parameters

        First are two examples of stopping tests.

        Note: task(1:4) must be assigned the value 'STOP' to terminate  
          the iteration and ensure that the final results are
          printed in the default format. The rest of the character
          string TASK may be used to store other information.

        1) Terminate if the total number of f and g evaluations
             exceeds MAXFEVAL.
      */
      if ( isave[33] >= MAXFEVAL )
        sprintf(task,"STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT") ;

      /*
        2) Terminate if  |proj g|/(1+|f|) < GRADSMALL, where 
           "proj g" denoted the projected gradient
      */
      if ( dsave[12] <= GRADSMALL * ( 1.0 + fabs(f) ) ) 
        sprintf(task,"STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL");

      /*
        We now wish to print the following information at each
        iteration (indexes correspond to fortran):
        
          1) the current iteration number, isave(30),
          2) the total number of f and g evaluations, isave(34),
          3) the value of the objective function f,
          4) the norm of the projected gradient,  dsve(13)

        See the comments at the end of driver1 for a description
        of the variables isave and dsave.
      */  
      
      printf("Iterate %d nfg = %d f = %e |proj g| = %e\n",isave[29],isave[33],f,
	     dsave[12]) ;      

      /*
        If the run is to be terminated, we print also the information
        contained in task as well as the final value of x.
      */
      if ( strncmp(task,"STOP",4) == 0 )
      {
        printf("%s\n",task) ;  
        printf("Final X = ") ;
        for ( i=0; i<n; i++ )
        {
          printf("  %11.4e",x[i]) ;
          if ( (i+1)%6 == 0 ) printf("\n") ;
        }
      }

      /* go back to the minimization routine. */
      lcont = 1 ;

    }

  } while ( lcont == 1 ) ;

  /* ---------- the end of the loop ------------- */
  free(iwa) ; iwa = NULL ;
  free(wa ) ; wa  = NULL ;
  free(g  ) ; g   = NULL ;
  free(xp ) ; xp  = NULL ;

  return (ans);
 
}

void gradfnum(int n, double *x, double *xp, double (* func)(int, double *, void*), 
              double *g, void *userdata)
{

  int i ;
  double fval, fvalp, perturb;
  
  for ( i = 0; i<n; i++ ) xp[i] = x[i] ;

  fval = func(n,x,userdata) ;
  for ( i = 0; i<n; i++ )
  {
    perturb = fabs(lbfgs_rt*x[i]);
    if ( perturb == 0.0 ) perturb = lbfgs_at ;
    xp[i] += perturb ;
    fvalp = func(n,xp,userdata) ;
    g[i] = (fvalp-fval)/perturb ;
    xp[i] -= perturb ;
  }

  return ;

}

void macheps()
{

  double sml = 1.0 ;
  while ( 1.0+sml != 1.0 ) sml *= 0.5 ;
  lbfgs_rt = sqrt(2.0*sml) ;
  lbfgs_at = lbfgs_rt ;

  return ;

}
