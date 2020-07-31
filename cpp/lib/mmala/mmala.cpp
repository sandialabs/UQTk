/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
// \file mala.cpp
// \author K. Sargsyan, C. Safta, B. Debusschere, 2012
// \brief MALA Markov chain Monte Carlo class

#include <math.h>
#include <float.h>
#include "error_handlers.h"
#include "deplapack.h"
#include "mala.h"
#include "mmala.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "mcmc.h"
#include "tmcmc.h"
#include "gen_defs.h"
#include "lbfgs_routines.h"

void MMALA::setMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *)){
  metricTensor_ = metricTensor;
  tensflag_ = true;
  return;
}

void MMALA::getMetricTensor(void (*metricTensor)(Array1D<double>&, Array2D<double>&, void *)){
  metricTensor = metricTensor_;
  return;
}

void MMALA::proposal(Array1D<double>& m_t,Array1D<double>& m_cand){
  int chol_info=0;
  char lu='L';

  Array1D<double> grads;
  this->gradlogPosterior_(m_t,grads,NULL);
  Array2D<double> mtensorinv;
  this->metricTensor_(m_t,mtensorinv,NULL);
  m_cand=m_t;

  Array1D<double> mtggrads;
  prodAlphaMatVec(mtensorinv, grads, 1.0, mtggrads) ;

  Array2D<double> sqrt_mtensorinv;
  sqrt_mtensorinv = mtensorinv;
  int chdim = this -> GetChainDim();
  FTN_NAME(dpotrf)(&lu,&chdim, sqrt_mtensorinv.GetArrayPointer(),&chdim,&chol_info);
  // Catch the error in Cholesky factorization
  if (chol_info != 0 )
    printf("Error in Cholesky factorization, info=%d\n", chol_info);

  for (int i=0; i < this -> GetChainDim(); i++) {
    m_cand(i) += this -> getEpsMALA() * this -> getEpsMALA() * mtggrads(i)/2.;
    for (int j=0; j < i+1; j++) {
      m_cand(i) += this -> getEpsMALA() *sqrt_mtensorinv(i,j)*dsfmt_genrand_nrv(&RandomState);
    }
  }

  return;
}

double MMALA::probOldNew(Array1D<double>& a, Array1D<double>& b){
  return 0.0;
  ///\todo In the original code there is a seperate branch for MMALA in the function probOldNew, but it is blank and has nothing in it. It would need to be defined in this setup.
}
