/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include "kl_utils.h"

#define XFILE "xgrid.dat"
#define YFILE "ygrid.dat"
#define GTYPE "cl0L"
#define AFAC   0.5
#define BFAC   1.1
#define NX     65
#define NY     65
#define DLEN   1.0

#define COVTYPE "SqExp"
#define CLEN    0.05
#define SIG     5.0

#define NKL   64
#define NSPL  128

/// \brief Displays information about this program
int usage(){

  printf("usage: kl_2D.x [-h] [-x] [-n<nx>] [-m<ny>] [-c<cov_type>] [-s<sigma>] [-l<cor_len>] [-e<nkl>] [-p<nspl>] [-w]\n");
  printf(" -h               : print out this help message \n");
  printf(" -x               : if used, read grids from xgrid.dat and ygrid.dat\n");
  printf(" -g <gtype>       : grid type (default=%s) \n",GTYPE);
  printf(" -a <afac>        : a-factor for double clustered grids (default=%e) \n",AFAC);
  printf(" -b <bfac>        : b-factor for clustered grids (default=%e) \n",BFAC);
  printf(" -n <nx>          : the number of x-grid points (default=%d) \n",NX);
  printf(" -m <ny>          : the number of y-grid points (default=%d) \n",NY);
  printf(" -c <cov_type>    : whether use analytical covariance (type: %s) or compute it from samples \n",COVTYPE);
  printf(" -s <sigma>       : standard deviation (default=%e) \n",SIG);
  printf(" -l <clen>        : correlation length (default=%e) \n",CLEN);
  printf(" -e <nkl>         : the number of KL modes retained (default=%d) \n",NKL);
  printf(" -p <nspl>        : the number of samples (default=%d) \n",NSPL);
  printf(" -w               : turn off saving samples to file \n");
  printf("================================================================================================================\n");
  printf("Input:: Nothing hard-coded.\n");
  printf("Output:: \n");
  printf("  - cov_out.dat:  covariance matrix\n");
  printf("  - eig.dat:      eigenvalues\n");
  printf("  - KLmodes.dat:  eigenmodes scaled with sqrt(eig)\n");
  printf("  - rel_diag.dat: pointwise covariance fraction explained by the finite KL expansion\n");
  printf("  - mean.dat:     pointwise mean based on samples\n");
  printf("  - xi_data.dat:  samples of eigenvalues\n");
  printf("================================================================================================================\n");
  exit(0);
  return 0;
}

/*!
   \brief Karhunen-Loeve decomposition of a UNIVARIATE process
   given covariance type or samples drawn from multivariate normal distribution
*/
int main(int argc, char *argv[])
{

  /* Set the default values */
  int    nkl  = NKL  ;
  int    nspl = NSPL ;

  int    nx   = NX   ;
  int    ny   = NY   ;
  double afac = AFAC ;
  double bfac = BFAC ;
  double dlen = DLEN ;

  bool  xflag = false;
  char* xfile = (char *) XFILE;
  char* yfile = (char *) YFILE;
  char *gtype = (char *) GTYPE;
  Array1D<double> xgrid, ygrid;

  bool   cflag = false ;
  double clen  = CLEN  ;
  double sigma = SIG   ;
  char* cov_type = (char *) COVTYPE;
  Array2D<double> cov ;

  bool   sflag = true ; /* by default save samples */

  /* Read the user input */
  int c;

  while ((c=getopt(argc,(char **)argv,"hxg:a:b:n:m:c:s:l:e:p:w"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'x':
      xflag = true;
      break;
    case 'g':
      gtype = optarg;
      break;
    case 'a':
      afac = strtod(optarg, (char **)NULL);
      break;
    case 'b':
      bfac = strtod(optarg, (char **)NULL);
      break;
    case 'n':
      nx   = strtol(optarg, (char **)NULL,0);
      break;
    case 'm':
      ny   = strtol(optarg, (char **)NULL,0);
      break;
    case 'c':
      cflag=true;
      cov_type = optarg;
      break;
    case 's':
      sigma = strtod(optarg, (char **)NULL);
      break;
   case 'l':
      clen = strtod(optarg, (char **)NULL);
      break;
    case 'e':
      nkl = strtol(optarg, (char **)NULL,0);
      break;
    case 'p':
      nspl = strtol(optarg, (char **)NULL,0);
      break;
    case 'w':
      sflag = false;
      break;
    }
  }

  int nxy=nx*ny;
  if ( nkl > nxy )
    throw Tantrum("kl_2D.cpp::main(): cannot request more KL modes than number of grid points");

  /* Print the input information on screen */
  cout << " - Number of grid points: " << nx <<"x"<<ny << endl<<flush;
  cout << " - Correlation length:    " << clen << endl<<flush;
  cout << " - Standard deviation:    " << sigma<< endl<<flush;
  cout << " - Number of KL modes:    " << nkl  << endl<<flush;
  if ( cflag )
    cout<<" - Will generate covariance of type "<<cov_type<<endl<<flush;
  else
    cout << " - Will generate covariance from "<<nspl<<" samples" << endl<<flush;
  if ( xflag )
    cout<<" - Will read grid from file: "<<xfile<<endl<<flush;
  else
    cout << " - Will create grid with "<<nx<<"x"<<ny<<" points" << endl<<flush;

  /* read grid from file or generate uniform grid on [0,1] */
  if ( xflag ) {
    /* read grid from file */
    read_datafile_1d(xgrid,xfile);
    read_datafile_1d(ygrid,yfile);
  }
  else {
    genGrid2D(xgrid,ygrid,nx,ny,dlen,gtype,afac,bfac);
    write_datafile_1d(xgrid,xfile);
    write_datafile_1d(ygrid,yfile);
  }
  Array2D<double> xygrid(nxy,2);
  for (int k=0; k < nxy; k++ ) {
    int i = k%nx;
    int j = k/nx;
    xygrid(k,0) = xgrid(i);
    xygrid(k,1) = ygrid(j);
  }

  Array2D<double> ySamples;

  cov.Resize(nxy,nxy,0.e0);
  if ( cflag ) {
    for ( int i = 0; i < nxy; i++)
      for ( int j = 0; j < nxy; j++)
        cov(i,j)=genCovAnl2D(xygrid,i,j,clen,sigma,cov_type);
  }
  // else {
  //   read_datafile(cov,"cov.dat");
  // }
  else {

    double dfac=1.0e-11;
    bool   tryagain = true;

    while ( ( tryagain ) && ( dfac < 1.e-6 ) ) {
      tryagain = false ;
      for ( int i = 0; i < nxy; i++)
        for ( int j = 0; j < nxy; j++)
          cov(i,j)=genCovAnl2D(xygrid,i,j,clen,sigma,string("SqExp"));
      for ( int i = 0; i < nxy; i++)
        cov(i,i) += dfac;
      //write_datafile( cov, "cov.dat" );

      /* Generate samples */
      char *lu = (char *) "L";
      int info ;
      FTN_NAME(dpotrf)( lu, &nxy, cov.GetArrayPointer(), &nxy, &info );

      /* Check the success in Cholesky factorization */
      if ( info != 0 ) {

        cout<<"Error in Cholesky factorization, info=" << info << endl << flush ;;
        dfac = dfac * 10.0 ;
        tryagain = true ;
        cout<<"  will try again with diagonal factor" << dfac << endl << flush ;;

      } /* done if Cholesky factorization fails */
    }
    dsfmt_t  rnstate ;
    int rseed=20120828;
    dsfmt_init_gen_rand(&rnstate, (uint32_t) rseed );

    Array1D<double> randSamples(nxy,0.0);
    ySamples.Resize(nxy,nspl,0.0);
    for ( int j = 0; j < nspl; j++) {
      for (int i = 0 ; i < nxy ; i++ )
        randSamples(i) = dsfmt_genrand_nrv(&rnstate);
      for ( int i = 0; i < nxy; i++ ) {
        ySamples(i,j)=0.0;
        for ( int k = 0; k < i+1; k++)
          ySamples(i,j) += (cov.GetArrayPointer())[i+k*nxy]*randSamples(k);
      }
    }
    if (sflag) write_datafile( ySamples, "samples.dat" );

    /* Compute samples mean */
    Array1D<double> mean(nxy,0.e0);
    for ( int i = 0 ; i < nxy ; i++ )
      mean(i) = getMean( ySamples, string("L"), i );
    write_datafile_1d( mean, "mean.dat" );

    /* Compute covariance matrix */

    /* 1. subtract the mean from samples */
    for ( int j = 0 ; j < nspl ; j++)
      for ( int i = 0 ; i < nxy ; i++)
        ySamples(i,j) -= mean(i) ;

    /* 2. compute the upper triangular part */
    for ( int i = 0; i < nxy; i++ ) {
      for ( int j = i; j < nxy; j++ ) {

        double dsum=0.0;
        for(int k = 0; k < nspl; k++ )
          dsum += ySamples(i,k)*ySamples(j,k);

        cov(i,j) = dsum/( (double) nspl );
      }
    }

    /* 3. transpose to fill out the lower triangle */
    for ( int i = 0; i < nxy; i++ )
      for ( int j = 0; j < i; j++ )
        cov(i,j) = cov(j,i) ;

  }
  write_datafile( cov, "cov.dat" );

  /*  Performing KL decomposition */
  cout << " --> Starting KL decomposition " << endl;

  /* Create 1D equivalent grid*/
  Array1D<double> xg1d(nxy,0.0);
  getGrid1dEquiv(xgrid,ygrid,xg1d);

  KLDecompUni decomposer(xg1d);
  int n_eig = decomposer.decompose(cov,nkl);

  if(n_eig <  nkl){
    printf("There are only %d  eigenvalues available (requested %d) \n",n_eig, nkl);
    nkl = n_eig;
  }

  const Array1D<double>& eigs    = decomposer.eigenvalues();
  const Array2D<double>& KLmodes = decomposer.KLmodes();

  cout << " --> KL decomposition done " << endl << flush;

  cout << "      - Obtained " << n_eig << " eigenvalues:" << endl;
  //for ( int i = 0; i < n_eig; i++)
  //  cout << "        " << i << " : " << eigs(i) << endl;
  write_datafile_1d(eigs,"eig.dat");

  cout << "      - Computing relative variances" << endl;
  Array1D<double> rel (nxy,0.0);
  for ( int i = 0; i < nxy; i++){
    for ( int j = 0; j < n_eig; j++){
      rel(i) += eigs(j)*KLmodes(i,j)*KLmodes(i,j);
    }
    rel(i) /= cov(i,i);
  }
  write_datafile_1d(rel,"relVar.dat");

  Array2D<double> scaledKLmodes(nxy,n_eig+2,0.0);
  for ( int k = 0; k < nxy; k++ ){
    int i = k%nx;
    int j = k/nx;
    scaledKLmodes(k,0) = xgrid(i);
    scaledKLmodes(k,1) = ygrid(j);
    for ( int l = 0; l < n_eig; l++ )
      scaledKLmodes(k,l+2) = KLmodes(k,l)*sqrt(eigs(l));
  }
  write_datafile(scaledKLmodes,"KLmodes.dat");

  if ( !cflag) {

    /* Project realizations onto KL modes */
    cout << "      - Project realizations onto KL modes " << endl;
    Array2D<double> xi(nkl, nspl, 0.e0);
    decomposer.KLproject( ySamples, xi );

    Array2D<double> xit(nspl,nkl,0.e0);
    transpose(xi,xit);
    write_datafile(xit,"xi_data.dat");

  }

  cout << " --> KL example done " << endl;

  return ( 0 ) ;

}
