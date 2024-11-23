/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
#include <cstdio>
#include <stddef.h>
#include <fstream>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "assert.h"
#include <getopt.h>

using namespace std;

#include "Array2D.h"

#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"
#include "error_handlers.h"
#include "ftndefs.h"
#include "depblas.h"
#include "deplapack.h"

double arrayMean(int n, double *a) ;
double arrayVar (int n, double *a, double aMean, int ddof) ;

#define VERB

void usage() {

  printf("Computes Sobol sensitivity indices via sampling\n");
  printf("usage: sens [-h] -d <ndim> -n<nspl> [-w<bw>] [-x<which_chaos>] [-a<alpha>] [-b<beta>]\n");
  printf(" -h          : print out this help message \n");
  printf(" -d <ndim>   : no. of dimensions (default=0) \n");
  printf(" -n <nspl>   : no. of samples (default=0) \n");
  printf(" -a <action> : action (default=NONE); can be set to spl* or idx*\n");
  printf("               where * can be FO (first order), TO (total order), or Jnt (joint)\n");
  printf(" -u <spl1>   : filename for first set of nspl samples (default=NONE)\n");
  printf(" -v <spl2>   : filename for second set of nspl samples (default=NONE)\n");
  printf(" -x <mev>    : filename for model evaluations (default=NONE)\n");
  printf(" -p <plist>  : filename for selected parameters (default=NONE)\n");
  exit(0);
  return;
}

int main(int argc, char *argv[]) {

  int ndim=0, nspl=0, incr=0 ;

  char *sfile   = (char *)"NONE" ;
  char *pfile   = (char *)"NONE" ; bool pflag =false ;
  char *m1file  = (char *)"NONE" ; bool m1flag=false ;
  char *m2file  = (char *)"NONE" ; bool m2flag=false ;
  char *action  = (char *)"NONE" ;

  dsfmt_t RandomState; /* random number structure */

  /* Read user input */
  int c;
  while ((c=getopt(argc,(char **)argv,"ha:x:p:u:v:d:n:"))!=-1){
    switch (c) {
    case 'h':
      usage();
      break;
    case 'a':
      action = optarg;
      break;
    case 'x':
      sfile = optarg;
      break;
    case 'p':
      pflag = true ;
      pfile = optarg ;
      break;
    case 'u':
      m1flag = true ;
      m1file = optarg ;
      break;
    case 'v':
      m2flag = true ;
      m2file = optarg ;
      break;
    case 'd':
      ndim = strtol(optarg, (char **)NULL,0); 
      break;
    case 'n':
      nspl = strtol(optarg, (char **)NULL,0); 
      break;
    }
  }

#ifdef VERB
  cout << "...No. of samples:    "<<nspl<<endl;
  cout << "...No. of dimensions: "<<ndim<<endl;
  cout << "...Action:            "<<action<<endl;
  cout << "...Sample filename:   "<<sfile<<endl;
  if ( pflag ) 
    cout << "Custom parameter associations filename: "<<pfile<<endl;
  if ( ( string(action) == string("splFO")  ) || 
       ( string(action) == string("splTO")  ) ||
       ( string(action) == string("splJnt") ) ) {
    cout << "...M1-matrix filename (if any): "<<m1file<<endl;
    cout << "...M2-matrix filename (if any): "<<m2file<<endl;
  }
#endif

  /* Initialize list of parameters */
  vector< vector<int> > plist;
  if ( pflag ){
    /* read list with custom parameter sub-lists */
    ifstream in( pfile );
    while ( in.good() ){
      int itmp;
      string theLine="";
      vector<int> pidx;
      getline(in,theLine);
      istringstream s(theLine);
      while( s >> itmp ) pidx.push_back(itmp); 
      if (pidx.size() > 0)
        plist.push_back(pidx);
    }
    in.close();
  } else {
    /* no custom list of parameters -> 
       all parameters are treated the same */
    for ( int j=0; j<ndim; j++ ) {
      vector<int> pidx(1,j);
      plist.push_back(pidx) ;
    }
  }

  /* Increment for scaling matrices */
  incr = 1;

  /* Read mat1 and mat2 matrices if necessary */
  Array2D<double> mat1, mat2, matj;
  double *pmat1=NULL, *pmat2=NULL;
  if ( ( string(action) == string("splFO")  ) || 
       ( string(action) == string("splTO")  ) ||
       ( string(action) == string("splJnt") ) ) {
    read_datafileVS(mat1,m1file);
    read_datafileVS(mat2,m2file);
    assert( ( (int)mat1.XSize() == nspl ) && ((int)mat1.YSize() == ndim) );
    assert( ( (int)mat2.XSize() == nspl ) && ((int)mat2.YSize() == ndim) );
    pmat1 = mat1.GetArrayPointer();
    pmat2 = mat2.GetArrayPointer();
  }

  if ( string(action) == string("splFO") ) {
       
    write_datafile(mat1, "splFO.txt");    
    matj = mat2;
    double *pmatj = matj.GetArrayPointer();
    for ( int j=0; j<plist.size(); j++ ) {
      for ( int i=0; i<plist[j].size(); i++ )
        FTN_NAME(dcopy)(&nspl,pmat1+plist[j][i]*nspl,&incr,pmatj+plist[j][i]*nspl,&incr);
      write_datafile(matj, "splFO.txt","a"); 
      for ( int i=0; i<plist[j].size(); i++ )
        FTN_NAME(dcopy)(&nspl,pmat2+plist[j][i]*nspl,&incr,pmatj+plist[j][i]*nspl,&incr);   
    }
    write_datafile(mat2, "splFO.txt","a");    

  }
 
  if ( string(action) == string("idxFO") ) {

    Array1D<double> ymat;
    read_datafileVS(ymat,sfile);
    int nidx = ymat.XSize()/nspl-2;

    Array1D<double> sFO(nidx);

    double *ymat1 = ymat.GetArrayPointer();
    double *ymat2 = ymat1+nspl*(nidx+1);

    /* Mean mat1*mat2 */
    Array1D<double> ym12(nspl,0.0);
    for ( int i=0; i<nspl; i++ ) ym12(i) = ymat1[i]*ymat2[i];
    double mean12 = arrayMean(nspl,ym12.GetArrayPointer());

    /* Mean & variance of mat1 */
    double mean1 = arrayMean( nspl, ymat1 );
    double vv1   = arrayVar ( nspl, ymat1, mean1, 1) ;

    for ( int j=1; j<nidx+1; j++ ) {
      double *ymatj = ymat1+j*nspl;
      for ( int i=0; i<nspl; i++ ) ym12(i) = ymat1[i]*ymatj[i];
      double vari = arrayMean ( nspl, ym12.GetArrayPointer()) ;
      sFO(j-1) = (vari*nspl/(nspl-1.0)-mean12)/vv1;
    }

    write_datafile_1d(sFO, "idxFO.txt");    

  }

  if ( string(action) == string("splTO") ) {
       
    write_datafile(mat1, "splTO.txt");    
    matj = mat1;
    double *pmatj = matj.GetArrayPointer();
    for ( int j=0; j<plist.size(); j++ ) {
      for ( int i=0; i<plist[j].size(); i++ )
        FTN_NAME(dcopy)(&nspl,pmat2+plist[j][i]*nspl,&incr,pmatj+plist[j][i]*nspl,&incr);
      write_datafile(matj, "splTO.txt","a"); 
      for ( int i=0; i<plist[j].size(); i++ )
        FTN_NAME(dcopy)(&nspl,pmat1+plist[j][i]*nspl,&incr,pmatj+plist[j][i]*nspl,&incr);
    }
    write_datafile(mat2, "splTO.txt","a");    

  }

  if ( string(action) == string("idxTO") ) {

    Array1D<double> ymat;
    read_datafileVS(ymat, sfile);
    int nidx = ymat.XSize()/nspl-2;

    Array1D<double> sTO(nidx);

    double *ymat1 = ymat.GetArrayPointer();
    double Ey  = arrayMean(nspl,ymat1);
    double vv1 = arrayVar (nspl,ymat1, Ey, 1) ;

    Array1D<double> ym12(nspl,0.0);
    for ( int j=1; j<nidx+1; j++ ) {

      double *ymatj = ymat1+j*nspl;
      for ( int i=0; i<nspl; i++ ) ym12(i) = pow(ymat1[i]-ymatj[i],2);
      sTO(j-1) = arrayMean ( nspl, ym12.GetArrayPointer())/(2.0*vv1);

    }

    write_datafile_1d(sTO, "idxTO.txt");    

  }

  if ( string(action) == string("splJnt") ) {
 
    /*
    for ( int k=0; k<plist.size(); k++ ) {
        cout<<"Line "<<k<<": Size "<<plist[k].size()<<":";
        for ( int i=0; i<plist[k].size(); i++ )
          cout<<plist[k][i]<<" ";
        cout<<endl;
    }
    */

    write_datafile(mat1, "splJ.txt");    
    matj = mat2;
    double *pmatj = matj.GetArrayPointer();
    for ( int k=0; k<plist.size()-1; k++ ) {
      for ( int j=k+1; j<plist.size(); j++ ) {
        for ( int i=0; i<plist[j].size(); i++ )
          FTN_NAME(dcopy)(&nspl,pmat1+plist[j][i]*nspl,&incr,pmatj+plist[j][i]*nspl,&incr);
        for ( int i=0; i<plist[k].size(); i++ )
          FTN_NAME(dcopy)(&nspl,pmat1+plist[k][i]*nspl,&incr,pmatj+plist[k][i]*nspl,&incr);
        write_datafile(matj, "splJ.txt", "a"); 
        for ( int i=0; i<plist[j].size(); i++ )
          FTN_NAME(dcopy)(&nspl,pmat2+plist[j][i]*nspl,&incr,pmatj+plist[j][i]*nspl,&incr);
        for ( int i=0; i<plist[k].size(); i++ )
          FTN_NAME(dcopy)(&nspl,pmat2+plist[k][i]*nspl,&incr,pmatj+plist[k][i]*nspl,&incr);
      }
    }
    write_datafile(mat2, "splJ.txt", "a");    

  }

  if ( string(action) == string("idxJnt") ) {

    Array1D<double> ymat;
    read_datafileVS(ymat,sfile);
    Array1D<double> sFO;
    read_datafileVS(sFO, (char *) "idxFO.txt");
    int nidx =  sFO.XSize();
    Array1D<double> sJnt(nidx*(nidx-1)/2);

    double *ymat1 = ymat.GetArrayPointer();
    double *ymat2 = ymat1+nspl*(nidx*(nidx-1)/2+1);

    /* Mean mat1*mat2 */
    Array1D<double> ym12(nspl,0.0);
    for ( int i=0; i<nspl; i++ ) ym12(i) = ymat1[i]*ymat2[i];
    double mean12 = arrayMean(nspl,ym12.GetArrayPointer());

    /* Mean & variance of mat1 */
    double mean1 = arrayMean( nspl, ymat1 );
    double vv1   = arrayVar ( nspl, ymat1, mean1, 1) ;

    int ijd=0;
    for ( int j=1; j<nidx+1; j++ ) {
      for ( int i=j+1; i<nidx+1; i++ ) {
        ijd++;
        double *ymatj = ymat1+ijd*nspl;
	      for ( int k=0; k<nspl; k++ ) ym12(k) = ymat1[k]*ymatj[k];
        double vari = arrayMean ( nspl, ym12.GetArrayPointer()) ;
        sJnt(ijd-1) = (vari*nspl/(nspl-1.0)-mean12)/vv1-sFO(j-1)-sFO(i-1);
      }
    }

    write_datafile_1d(sJnt, "idxJnt.txt");    

  }

  return (0) ;

}

double arrayMean(int n, double *a) {

  int    i, m    ;
  double dsum = 0.0 ;
  
  m = n % 5;
  if ( m != 0 ) 
    for ( i=0; i < m ; i++ ) dsum += a[i];

  if ( n<5 ) return (dsum/n);

  for ( i=m; i < n ; i+=5 ) {
    dsum += a[i]+a[i+1]+a[i+2]+a[i+3]+a[i+4];
  }
  return (dsum/n);
  
}

double arrayVar(int n, double *a, double aMean, int ddof) {

  int    i, m    ;
  double dsum = 0.0 ;
  
  m = n % 5;
  if ( m != 0 ) 
    for ( i=0; i < m ; i++ ) dsum +=  pow(a[i]-aMean,2);

  if ( n<5 ) return (dsum/(n-ddof));

  for ( i=m; i < n ; i+=5 ) {
    dsum += pow(a[i  ]-aMean,2); 
    dsum += pow(a[i+1]-aMean,2);
    dsum += pow(a[i+2]-aMean,2);
    dsum += pow(a[i+3]-aMean,2);
    dsum += pow(a[i+4]-aMean,2);
  }
  return (dsum/(n-ddof));
  
}
