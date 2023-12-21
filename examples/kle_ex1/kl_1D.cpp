/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.4
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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#include "kl_utils.h"

#define XFILE   "xgrid.dat"
#define GTYPE   "cl0"
#define AFAC    0.5
#define BFAC    1.1
#define NPTS    129
#define DLEN    1.0

#define COVTYPE "SqExp"
#define CLEN  0.05
#define SIG   5.0

#define NKL   64
#define NSPL  128


/// \brief Displays information about this program
int usage(){

  printf("usage: kl_1D.x [-h]  [-x<grid_file>] [-g<gtype>] [-a<afac>] [-b<bfac>] [-n<npts>] [-c<cov_type>] [-s<sigma>] [-l<cor_len>][-e<nkl>] [-p<nspl>] [-w]\n");
  printf(" -h               : print out this help message \n");
  printf(" -x <xfile>       : define the file from which the time grid is being read (default=%s) \n",XFILE);
  printf(" -g <gtype>       : grid type (default=%s) \n",GTYPE);
  printf(" -a <afac>        : a-factor for double clustered grids (default=%e) \n",AFAC);
  printf(" -b <bfac>        : b-factor for clustered grids (default=%e) \n",BFAC);
  printf(" -n <npts>        : the number of grid points (default=%d) \n",NPTS);
  printf(" -c <cov_type>    : wether use analytical covariance (default type: %s) or compute it from samples \n",COVTYPE);
  printf(" -s <sigma>       : standard deviation (default=%e) \n",SIG);
  printf(" -l <clen>        : correlation length (default=%e) \n",CLEN);
  printf(" -e <nkl>         : the number of KL modes retained (default=%d) \n",NKL);
  printf(" -p <nspl>        : the number of samples (default=%d) \n",NSPL);
  printf(" -w               : turn off saving samples to file \n");
  printf("================================================================================================================\n");
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

  int    npts = NPTS ;
  double afac = AFAC ;
  double bfac = BFAC ;
  double dlen = DLEN ;
  bool  xflag = false;
  char *xfile = (char *) XFILE;
  char *gtype = (char *) GTYPE;
  Array1D<double> xgrid;

  bool   cflag = false ;
  double clen  = CLEN  ;
  double sigma = SIG   ;
  char* cov_type = (char *) COVTYPE;
  Array2D<double> cov ;

  bool   sflag = true ; /* by default save samples */

  /* Read the user input */
  int c;

  while ((c=getopt(argc,(char **)argv,"hx:g:a:b:n:d:c:s:l:e:p:w"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'x':
      xflag = true;
      xfile = optarg;
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
      npts = strtol(optarg, (char **)NULL,0);
      break;
    case 'd':
      dlen = strtod(optarg, (char **)NULL);
      break;
    case 'c':
      cflag = true;
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

  if ( nkl > npts )
    throw Tantrum("kl_1D.cpp::main(): cannot request more KL modes than number of grid points");

  /* Print the input information on screen */
  cout << " - Number of grid points: " << npts << endl<<flush;
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
    cout << " - Will create grid with "<<npts<<" points in [0,1]" << endl<<flush;

  /* read grid from file or generate uniform grid on [0,1] */
  if ( xflag )
    /* read grid from file */
    read_datafile_1d(xgrid,xfile);
  else {
    genGrid1D(xgrid,npts,dlen,gtype,afac,bfac);
    write_datafile_1d(xgrid,xfile);
  }
  Array2D<double> ySamples;

  cov.Resize(npts,npts,0.e0);
  if ( cflag ) {
    for ( int i = 0; i < npts; i++)
      for ( int j = 0; j < npts; j++)
        cov(i,j)=genCovAnl1D(xgrid(i),xgrid(j),clen,sigma,cov_type);
  }
  else {

    for ( int i = 0; i < npts; i++)
      for ( int j = 0; j < npts; j++)
        cov(i,j)=genCovAnl1D(xgrid(i),xgrid(j),clen,sigma,string("SqExp"));
    for ( int i = 0; i < npts; i++)
      cov(i,i) += 1.e-13;

    /* Generate samples */
    char *lu = (char *) "L";
    int info ;
    FTN_NAME(dpotrf)( lu, &npts, cov.GetArrayPointer(), &npts, &info );

    /* Check the success in Cholesky factorization */
    if ( info != 0 ) {

      cout<<"Error in Cholesky factorization, info=" << info << endl << flush ;;
      exit(1);

    } /* done if Cholesky factorization fails */

    dsfmt_t  rnstate ;
    int rseed = 20120828;
    dsfmt_init_gen_rand(&rnstate, (uint32_t) rseed );

    float progress = 0.0;
    int barWidth = 70;
    cout << " - Generate samples" << endl<<flush;
    Array1D<double> randSamples(npts,0.0);
    ySamples.Resize(npts,nspl,0.0);
    for ( int j = 0; j < nspl; j++) {

      if (float(j) / nspl > progress+0.01) {

        cout << "   [";
        int pos = barWidth * progress;
        for (int ii = 0; ii< barWidth; ++ii) {
            if (ii < pos) cout << "=";
            else if (ii == pos) cout << ">";
            else cout << " ";
        }
        cout << "] " << int(progress * 100.0) << " %\r";
        cout << flush;
        progress += 0.01;
      }

      for (int i = 0 ; i < npts ; i++ )
        randSamples(i) = dsfmt_genrand_nrv(&rnstate);
      for ( int i = 0; i < npts; i++ ) {
        ySamples(i,j)=0.0;
        for ( int k = 0; k < i+1; k++)
          ySamples(i,j) += (cov.GetArrayPointer())[i+k*npts]*randSamples(k);
      }
    }

    if (sflag) write_datafile( ySamples, "samples.dat" );

    /* Compute samples mean */
    Array1D<double> mean(npts,0.e0);
    string dirtype("L");
    for ( int i = 0 ; i < npts ; i++ )
      mean(i) = getMean( ySamples, dirtype, i );
    write_datafile_1d( mean, "mean.dat" );

    /* Compute covariance matrix */

    /* 1. subtract the mean from samples */
    for ( int j = 0 ; j < nspl ; j++)
      for ( int i = 0 ; i < npts ; i++)
        ySamples(i,j) -= mean(i) ;

    /* 2. compute the upper triangular part */
    for ( int i = 0; i < npts; i++ ) {
      for ( int j = i; j < npts; j++ ) {

        double dsum=0.0;
        for(int k = 0; k < nspl; k++ )
          dsum += ySamples(i,k)*ySamples(j,k);

        cov(i,j) = dsum/( (double) nspl );
      }
    }

    /* 3. transpose to fill out the lower triangle */
    for ( int i = 0; i < npts; i++ )
      for ( int j = 0; j < i; j++ )
        cov(i,j) = cov(j,i) ;

  }
  write_datafile( cov, "cov.dat" );

  /*  Performing KL decomposition */
  cout << " --> Starting KL decomposition " << endl;

  KLDecompUni decomposer(xgrid);
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
  Array1D<double> rel (npts,0.0);
  for ( int i = 0; i < npts; i++){
    for ( int j = 0; j < n_eig; j++){
      rel(i) += eigs(j)*KLmodes(i,j)*KLmodes(i,j);
    }
    rel(i) /= cov(i,i);
  }
  write_datafile_1d(rel,"relVar.dat");

  Array2D<double> scaledKLmodes(npts,n_eig+1,0.0);
  for ( int i = 0; i < npts; i++ ){
    scaledKLmodes(i,0) = xgrid(i);
    for ( int j = 0; j < n_eig; j++ )
      scaledKLmodes(i,j+1) = KLmodes(i,j)*sqrt(eigs(j));
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
