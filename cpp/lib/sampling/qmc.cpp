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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#include "sampling.hpp"

/**
 \brief Returns first nPrime prime numbers
 */
int Sampling::getPrimes(const int nPrime, std::vector<int> &primes) {

  if ( nPrime < 1 ) return (-1) ; // return error if no. of primes is less than 1

  primes[0] = 2;

  int pTest = 3;
  for ( int i = 2 ; i <=nPrime ; ) {

    int status = 1;
    for ( int divN = 2 ; divN <= (int)sqrt(pTest) ; divN++ ) {
      if ( pTest%divN == 0 ) {
        status = 0;
        break;
      }
    }
    if ( status != 0 ) {
      primes[i-1] = pTest;
      i++;
    }
    status = 1;
    pTest++;
   }

  return (0);

}

/**
 \brief Returns nelem elements, dim-dimensional, Halton QMC sequence
 */
void Sampling::getHaltonSeq ( const int nelem, const int dim, std::vector<double> &seq ) {


  this->getHaltonSeq( nelem, dim, this->hStep_, this->hSkip_, this->hJump_, this->hBase_, seq );

  this->hStep_ += nelem;

  return;

}

/**
 \brief Returns nelem elements, dim-dimensional, Halton QMC sequence
 */
void Sampling::getHaltonSeq ( const int nelem, const int dim, const int step,
		    const std::vector<int> &skip, const std::vector<int> &jump, const std::vector<int> &base,
		    std::vector<double> &seq) {

  assert(nelem*dim<=seq.size());

  for ( int i=0; i<nelem*dim; i++ ) seq[0]=0.0;

  for ( int i = 0; i < dim; i++ ) {
    for ( int j = 0; j < nelem; j++ ) {

      int seq1D = skip[i] + ( step + j ) * jump[i];
      double obase = 1.0 / ( ( double ) base[i] );

      while ( seq1D != 0 ) {
        int digit = seq1D % base[i];
        seq[j*dim+i] += ( ( double ) digit ) * obase ;
        obase = obase / ( ( double ) base[i] );
        seq1D = seq1D / base[i];
      }
    }
  }

  return;

}

/**
 \brief Returns nelem elements, dim-dimensional, Hammersley QMC sequence
 */
void Sampling::getHammersleySeq ( const int nelem, const int dim, double *seq ) {

  for ( int i=0; i<nelem*dim; i++ ) seq[0]=0.0;

  for ( int j = 0; j < nelem; j++ )
    seq[j*dim] = ( double ) ( j % ( nelem + 1 ) ) / ( double ) ( nelem );

  for ( int i = 1; i < dim; i++ ) {

    for ( int j = 0; j < nelem; j++ ) {

      int seq1D = j;
      double obase = 1.0 / ( ( double ) this->hBase_[i-1] );

      while ( seq1D != 0 ) {
        int digit = seq1D % this->hBase_[i-1];
        seq[j*dim+i] += ( ( double ) digit ) * obase ;
        obase = obase / ( ( double ) this->hBase_[i-1] );
        seq1D = seq1D / this->hBase_[i-1];
      }
    }
  }

  return ;

}
