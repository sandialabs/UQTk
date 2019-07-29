#include "dsfmt_add.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"


#define DPI2 8.0*atan(1.0)

static int dsfmt_nrv_calls = 0 ;

void dsfmt_reset_add()
{
  dsfmt_nrv_calls = 0 ;
}

double dsfmt_gv_genrand_urv()
{
 return dsfmt_gv_genrand_open_open();
}

double dsfmt_genrand_urv(dsfmt_t *dsfmt)
{
 return dsfmt_genrand_open_open(dsfmt);
}


double dsfmt_gv_genrand_urv_sm(double a, double b)
{
  if (a>b){
    printf("dsfmt_gv_genrand_unif_sm(): Lower bound of the interval is larger than the upper bound : %lg > %lg\n",a, b);
    printf(" -> Abort !\n");
    fflush(stdout);
    exit(0) ;
  }

  return ( dsfmt_gv_genrand_urv()*(b-a)+a );
}

double dsfmt_genrand_urv_sm(dsfmt_t *dsfmt, double a, double b)
{
 if (a>b){
    printf("dsfmt_genrand_unif_sm(): Lower bound of the interval is larger than the upper bound : %lg > %lg\n",a, b);
    printf(" -> Abort !\n");
    fflush(stdout);
    exit(0) ;
  }

  return ( dsfmt_genrand_urv(dsfmt)*(b-a)+a );
}

double dsfmt_gv_genrand_nrv()
{
  double unif[2], sm2logu1, pi2u2 ;
  static double zn1, zn2;
  if (  dsfmt_nrv_calls == 0 )
  {
    unif[0]=dsfmt_gv_genrand_open_close();
    unif[1]=dsfmt_gv_genrand_open_close();
    sm2logu1 = sqrt(-2.0*log(unif[0])) ;
    pi2u2    = DPI2*unif[1] ;
    zn1 = sm2logu1*cos(pi2u2) ;
    zn2 = sm2logu1*sin(pi2u2) ;
    dsfmt_nrv_calls = 1;
    return ( zn1 ) ;
  }
  else if (  dsfmt_nrv_calls == 1 )
  {
    dsfmt_nrv_calls = 0 ;
    return ( zn2 );
  }
  else
  {
    printf("dsfmt_genrand_nrv(): Unknown value for internal variable dsfmt_nrv_calls : %d\n",dsfmt_nrv_calls);
    printf(" -> Abort !\n");
    fflush(stdout);
    exit(0) ;
  }

}

double dsfmt_genrand_nrv(dsfmt_t *dsfmt)
{
  double unif[2], sm2logu1, pi2u2 ;
  static double zn1, zn2;

  // if (  dsfmt_nrv_calls == 0 )
  // {
    unif[0] = dsfmt_genrand_open_close(dsfmt);
    unif[1] = dsfmt_genrand_open_close(dsfmt);
    sm2logu1 = sqrt(-2.0*log(unif[0])) ;
    pi2u2    = DPI2*unif[1] ;
    zn1 = sm2logu1*cos(pi2u2) ;
    // zn2 = sm2logu1*sin(pi2u2) ;
    // dsfmt_nrv_calls = 1;
    return ( zn1 ) ;
  // }
  // else if (  dsfmt_nrv_calls == 1 )
  // {
  //   dsfmt_nrv_calls = 0 ;
  //   return ( zn2 );
  // }
  // else
  // {
  //   printf("dsfmt_genrand_nrv(): Unknown value for internal variable dsfmt_nrv_calls : %d\n",dsfmt_nrv_calls);
  //   printf(" -> Abort !\n");
  //   fflush(stdout);
  //   exit(0) ;
  // }

}

double dsfmt_gv_genrand_nrv_sm(double mu, double sigma)
{
  if (sigma <= 0) {
    printf("dsfmt_gv_genrand_nrv_sm(): Sigma is less than or equal to 0");
    printf(" -> Abort !\n");
    fflush(stdout);
    exit(0) ;
  }
  return ( dsfmt_gv_genrand_nrv()*sigma+mu );
}

double dsfmt_genrand_nrv_sm(dsfmt_t *dsfmt, double mu, double sigma)
{
  if (sigma <= 0) {
    printf("dsfmt_genrand_nrv(): Sigma is less than or equal to 0");
    printf(" -> Abort !\n");
    fflush(stdout);
    exit(0);
  }
  return ( dsfmt_genrand_nrv(dsfmt)*sigma+mu );
}

