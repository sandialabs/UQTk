#ifndef DSFMTADD_H_Seen
#define DSFMTADD_H_Seen

#include "dSFMT.h"

#ifdef __cplusplus
extern "C"
{
#endif

  void dsfmt_reset_add() ;

  double dsfmt_gv_genrand_urv();
  double dsfmt_genrand_urv(dsfmt_t *dsfmt);
  double dsfmt_gv_genrand_urv_sm(double a, double b);
  double dsfmt_genrand_urv_sm(dsfmt_t *dsfmt, double a, double b);
  
  double dsfmt_gv_genrand_nrv();
  double dsfmt_genrand_nrv(dsfmt_t *dsfmt) ;
  double dsfmt_gv_genrand_nrv_sm(double mu, double sigma) ;
  double dsfmt_genrand_nrv_sm(dsfmt_t *dsfmt, double mu, double sigma) ;
#ifdef __cplusplus
}
#endif

#endif
