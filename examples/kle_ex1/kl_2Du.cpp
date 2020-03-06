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
#include "kl_utils.h"

#define NKL   64
#define NSPL  128
#define CLEN  0.05
#define SIG   5.0

#define COVTYPE "SqExp"
#define XFILE   "data/cali_grid.dat"
#define TFILE   "data/cali_tria.dat"


/// \brief Displays information about this program
int usage(){

  printf("usage: cor_kl [-h]  [-c<cov_type>] [-e<nkl>] [-p<nspl>] [-s<sigma>] [-l<cor_len>] [-x<xfile>] [-t<tfile>] [-w]\n");
  printf(" -h               : print out this help message \n");
  printf(" -c <cov_type>    : whether use analytical covariance (type: %s) or compute it from samples \n",COVTYPE);
  printf(" -e <nkl>         : number of KL modes retained (default=%d) \n",NKL);
  printf(" -p <nspl>        : number of samples (default=%d) \n",NSPL);
  printf(" -s <sigma>       : standard deviation (default=%e) \n",SIG);
  printf(" -l <clen>        : correlation length (default=%e) \n",CLEN);
  printf(" -w               : turn off saving samples to file \n");
  printf("================================================================================================\n");
  printf("Output:: \n");
  printf("  - cov_out.dat:  covariance matrix\n");
  printf("  - eig.dat:      eigenvalues\n");
  printf("  - KLmodes.dat:  eigenmodes scaled with sqrt(eig)\n");
  printf("  - rel_diag.dat: pointwise covariance fraction explained by the finite KL expansion\n");
  printf("  - mean.dat:     pointwise mean based on samples\n");
  printf("  - xi_data.dat:  samples of eigenvalues\n");
  printf("================================================================================================\n");
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

  char *xfile = (char *) XFILE;
  char *tfile = (char *) TFILE;
  Array2D<double> xgrid;
  Array2D<int>    tgrid;

  bool   cflag = false ;
  double clen  = CLEN  ;
  double sigma = SIG   ;
  char *cov_type = (char *) COVTYPE;
  Array2D<double> cov ;

  bool   sflag = true ; /* by default save samples */

  /* Read the user input */
  int c;

  while ((c=getopt(argc,(char **)argv,"hc:e:p:s:l:t:x:w"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'c':
      cflag=true;
      cov_type = optarg;
      break;
    case 'e':
      nkl = strtol(optarg, (char **)NULL,0);
      break;
    case 'p':
      nspl = strtol(optarg, (char **)NULL,0);
      break;
    case 's':
      sigma = strtod(optarg, (char **)NULL);
      break;
   case 'l':
      clen = strtod(optarg, (char **)NULL);
      break;
    case 'x':
      xfile = optarg;
      break;
    case 't':
      tfile = optarg;
      break;
    case 'w':
      sflag = false;
      break;
    }
  }

  /* Print the input information on screen */
  cout << " - Correlation length:    " << clen << endl<<flush;
  cout << " - Standard deviation:    " << sigma<< endl<<flush;
  cout << " - Number of KL modes:    " << nkl  << endl<<flush;
  if ( cflag )
    cout<<" - Will generate covariance of type "<<cov_type<<endl<<flush;
  else
    cout << " - Will generate covariance from "<<nspl<<" samples" << endl<<flush;

  /* read grid from file */
  read_datafileVS(xgrid,xfile);
  read_datafileVS(tgrid,tfile);
  cout<<" - No. of grid points:    " <<xgrid.XSize() << endl << flush ;
  cout<<" - No. of triangles  :    " <<tgrid.XSize() << endl << flush ;

  int nxy = (int) xgrid.XSize() ;
  if ( nkl > nxy )
    throw Tantrum("kl_sample2Du::main(): cannot request more KL modes than number of grid points");

  Array2D<double> ySamples;

  cov.Resize(nxy,nxy,0.e0);
  if ( cflag ) {
    for ( int i = 0; i < nxy; i++)
      for ( int j = 0; j < nxy; j++)
        cov(i,j)=genCovAnl2D(xgrid,i,j,clen,sigma,cov_type);
  }
  else {

    double dfac=1.0e-12;
    bool   tryagain = true;

    while ( ( tryagain ) && ( dfac < 1.e-6 ) ) {

      tryagain = false ;
      for ( int i = 0; i < nxy; i++)
        for ( int j = 0; j < nxy; j++)
          cov(i,j)=genCovAnl2D(xgrid,i,j,clen,sigma,string("SqExp"));
      for ( int i = 0; i < nxy; i++) cov(i,i) += dfac;
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
  Array1D<double> weights(nxy,0.0);
  getWeights2DU(xgrid,tgrid,weights);

  KLDecompUni decomposer;
  decomposer.SetWeights(weights);
  int n_eig = decomposer.decompose(cov,nkl);

  if(n_eig <  nkl){
    printf("There are only %d  eigenvalues available (requested %d) \n",n_eig, nkl);
    nkl = n_eig;
  }

  const Array1D<double>& eigs    = decomposer.eigenvalues();
  const Array2D<double>& KLmodes = decomposer.KLmodes();

  cout << " --> KL decomposition done " << endl << flush;

  cout << "      - Obtained " << n_eig << " eigenvalues:" << endl;
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
    scaledKLmodes(k,0) = xgrid(k,0);
    scaledKLmodes(k,1) = xgrid(k,1);
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
