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
#include "math.h"
#include "tools.h"
#include "gkplib.h"

#define MAX(a,b)        (((a) > (b)) ? (a) : (b))


/* (1) */
static double x1[] = {0.0000000};
static double w1[] = {2.0000000};

/* (1+2) */
static double x3[] = {-0.77459666924148337704,0.0,                    0.77459666924148337704 };
static double w3[] = {0.555555555555555555556,0.888888888888888888889,0.555555555555555555556};

/* (1+2+4) */
static double x7[] = {-0.96049126870802028342,-0.77459666924148337704,-0.43424374934680255800,
                       0.0,
                       0.43424374934680255800, 0.77459666924148337704, 0.96049126870802028342};
static double w7[] = { 0.104656226026467265194,0.268488089868333440729,0.401397414775962222905,
                       0.450916538658474142345,
                       0.401397414775962222905,0.268488089868333440729,0.104656226026467265194};

/* (1+2+4+8) */
static double x15[] = {-0.99383196321275502221,-0.96049126870802028342,-0.88845923287225699889,
                       -0.77459666924148337704,-0.62110294673722640294,-0.43424374934680255800,
                       -0.22338668642896688163, 0.0,                    0.22338668642896688163,
                        0.43424374934680255800, 0.62110294673722640294, 0.77459666924148337704,
                        0.88845923287225699889, 0.96049126870802028342, 0.99383196321275502221 };
static double w15[] = {0.0170017196299402603390,0.0516032829970797396969,0.0929271953151245376859,
                       0.134415255243784220360, 0.171511909136391380787, 0.200628529376989021034,
                       0.219156858401587496404, 0.225510499798206687386, 0.219156858401587496404,
                       0.200628529376989021034, 0.171511909136391380787, 0.134415255243784220360,
                       0.0929271953151245376859,0.0516032829970797396969, 0.0170017196299402603390};

/* (1+2+4+8+16) */
static double x31[] = {-0.99909812496766759766,-0.99383196321275502221,-0.98153114955374010687,
                       -0.96049126870802028342,-0.92965485742974005667,-0.88845923287225699889,
                       -0.83672593816886873550,-0.77459666924148337704,-0.70249620649152707861,
                       -0.62110294673722640294,-0.53131974364437562397,-0.43424374934680255800,
                       -0.33113539325797683309,-0.22338668642896688163,-0.11248894313318662575,
                        0.0,
                        0.11248894313318662575, 0.22338668642896688163, 0.33113539325797683309,
                        0.43424374934680255800, 0.53131974364437562397, 0.62110294673722640294,
                        0.70249620649152707861, 0.77459666924148337704, 0.83672593816886873550,
                        0.88845923287225699889, 0.92965485742974005667, 0.96049126870802028342,
                        0.98153114955374010687, 0.99383196321275502221, 0.99909812496766759766 };
static double w31[] = {0.00254478079156187441540,0.00843456573932110624631,0.0164460498543878109338,
                       0.0258075980961766535646, 0.0359571033071293220968, 0.0464628932617579865414,
                       0.0569795094941233574122, 0.0672077542959907035404, 0.0768796204990035310427,
                       0.0857559200499903511542, 0.0936271099812644736167, 0.100314278611795578771,
                       0.105669893580234809744,  0.109578421055924638237,  0.111956873020953456880,
                       0.112755256720768691607,
                       0.111956873020953456880,  0.109578421055924638237,  0.105669893580234809744,
                       0.100314278611795578771,  0.0936271099812644736167, 0.0857559200499903511542,
                       0.0768796204990035310427, 0.0672077542959907035404, 0.0569795094941233574122,
                       0.0464628932617579865414, 0.0359571033071293220968, 0.0258075980961766535646,
                       0.0164460498543878109338, 0.00843456573932110624631,0.00254478079156187441540 };

/* (1+2+4+8+16+32) */
static double x63[] = {-0.99987288812035761194,-0.99909812496766759766,-0.99720625937222195908,
                       -0.99383196321275502221,-0.98868475754742947994,-0.98153114955374010687,
                       -0.97218287474858179658,-0.96049126870802028342,-0.94634285837340290515,
                       -0.92965485742974005667,-0.91037115695700429250,-0.88845923287225699889,
                       -0.86390793819369047715,-0.83672593816886873550,-0.80694053195021761186,
                       -0.77459666924148337704,-0.73975604435269475868,-0.70249620649152707861,
                       -0.66290966002478059546,-0.62110294673722640294,-0.57719571005204581484,
                       -0.53131974364437562397,-0.48361802694584102756,-0.43424374934680255800,
                       -0.38335932419873034692,-0.33113539325797683309,-0.27774982202182431507,
                       -0.22338668642896688163,-0.16823525155220746498,-0.11248894313318662575,
                       -0.056344313046592789972,0.0,                    0.056344313046592789972,
                        0.11248894313318662575, 0.16823525155220746498, 0.22338668642896688163,
                        0.27774982202182431507, 0.33113539325797683309, 0.38335932419873034692,
                        0.43424374934680255800, 0.48361802694584102756, 0.53131974364437562397,
                        0.57719571005204581484, 0.62110294673722640294, 0.66290966002478059546,
                        0.70249620649152707861, 0.73975604435269475868, 0.77459666924148337704,
                        0.80694053195021761186, 0.83672593816886873550, 0.86390793819369047715,
                        0.88845923287225699889, 0.91037115695700429250, 0.92965485742974005667,
                        0.94634285837340290515, 0.96049126870802028342, 0.97218287474858179658,
                        0.98153114955374010687, 0.98868475754742947994, 0.99383196321275502221,
                        0.99720625937222195908, 0.99909812496766759766, 0.99987288812035761194 };
static double w63[] = {0.000363221481845530659694,0.00126515655623006801137,0.00257904979468568827243,
                       0.00421763044155885483908, 0.00611550682211724633968,0.00822300795723592966926,
                       0.0104982469096213218983,  0.0129038001003512656260, 0.0154067504665594978021,
                       0.0179785515681282703329,  0.0205942339159127111492, 0.0232314466399102694433,
                       0.0258696793272147469108,  0.0284897547458335486125, 0.0310735511116879648799,
                       0.0336038771482077305417,  0.0360644327807825726401, 0.0384398102494555320386,
                       0.0407155101169443189339,  0.0428779600250077344929, 0.0449145316536321974143,
                       0.0468135549906280124026,  0.0485643304066731987159, 0.0501571393058995374137,
                       0.0515832539520484587768,  0.0528349467901165198621, 0.0539054993352660639269,
                       0.0547892105279628650322,  0.0554814043565593639878, 0.0559784365104763194076,
                       0.0562776998312543012726,  0.0563776283603847173877, 0.0562776998312543012726,  
                       0.0559784365104763194076,  0.0554814043565593639878, 0.0547892105279628650322,
                       0.0539054993352660639269,  0.0528349467901165198621, 0.0515832539520484587768,
                       0.0501571393058995374137,  0.0485643304066731987159, 0.0468135549906280124026,
                       0.0449145316536321974143,  0.0428779600250077344929, 0.0407155101169443189339,
                       0.0384398102494555320386,  0.0360644327807825726401, 0.0336038771482077305417,
                       0.0310735511116879648799,  0.0284897547458335486125, 0.0258696793272147469108,
                       0.0232314466399102694433,  0.0205942339159127111492, 0.0179785515681282703329,
                       0.0154067504665594978021,  0.0129038001003512656260, 0.0104982469096213218983,
                       0.00822300795723592966926, 0.00611550682211724633968,0.00421763044155885483908,
                       0.00257904979468568827243, 0.00126515655623006801137,0.000363221481845530659694 };

/* (1) */
static double xn1[] = {0.0000000000000000};
static double wn1[] = {1.0000000000000000};

/* (1+2) */
static double xn3[] = {-1.73205080756887719, 0.000000000000000000, 1.73205080756887719};
static double wn3[] = {0.166666666666666657, 0.66666666666666663, 0.166666666666666657};

/* (1+2+6) */
static double xn9[] = {-4.18495601767273229, -2.86127957605705818, -1.73205080756887719, 
                       -0.741095349994540853, 0.00000000000000000,  0.741095349994540853, 
                        1.73205080756887719,  2.86127957605705818,  4.18495601767273229 };
static double wn9[] = { 9.42694575565174701E-05, 0.00799632547089352934, 0.0948509485094851251,  
			0.270074329577937755, 0.253968253968254065, 0.270074329577937755,  
                        0.0948509485094851251,0.00799632547089352934,9.42694575565174701E-05 };

/* (1+2+6+10) */
static double xn19[] = {-6.36339449433636961,  -5.18701603991365623, -4.18495601767273229,  
                        -3.20533379449919442,  -2.86127957605705818, -2.59608311504920231,  
                        -1.73205080756887719,  -1.23042363402730603, -0.741095349994540853,  
                         0.0000000000000000,      
                         0.741095349994540853,  1.23042363402730603,  1.73205080756887719,     
			 2.59608311504920231,    2.86127957605705818,  3.20533379449919442,     
			 4.18495601767273229,    5.18701603991365623,  6.36339449433636961 };
static double wn19[] = { 8.62968460222986318E-10, 6.09480873146898402E-07, 6.01233694598479965E-05, 
                         0.00288488043650675591, -0.00633722479337375712,  0.0180852342547984622,  
                         0.0640960546868076103,   0.0611517301252477163,   0.208324991649608771,  
                         0.303467199854206227,    
                         0.208324991649608771,    0.0611517301252477163,  0.0640960546868076103, 
			 0.0180852342547984622,  -0.00633722479337375712, 0.00288488043650675591,  
			 6.01233694598479965E-05, 6.09480873146898402E-07,8.62968460222986318E-10 };

/* (1+2+6+10+16) */
static double xn35[] = {-9.0169397898903032,   -7.98077179859056063,  -7.12210670080461661, 
                        -6.36339449433636961,  -5.69817776848810986,  -5.18701603991365623, 
                        -4.73643308595229673,  -4.18495601767273229,  -3.63531851903727832, 
                        -3.20533379449919442,  -2.86127957605705818,  -2.59608311504920231, 
                        -2.23362606167694189,  -1.73205080756887719,  -1.23042363402730603, 
                        -0.741095349994540853, -0.248992297579960609,  
                         0.00000000000000000, 
                         0.248992297579960609,  0.741095349994540853,  
                         1.23042363402730603,   1.73205080756887719,   2.23362606167694189, 
                         2.59608311504920231,   2.86127957605705818,   3.20533379449919442,    
                         3.63531851903727832,   4.18495601767273229,   4.73643308595229673,    
                         5.18701603991365623,   5.69817776848810986,   6.36339449433636961,     
                         7.12210670080461661,   7.98077179859056063,   9.0169397898903032 };
static double wn35[] = { 1.05413265823340136E-18, 5.45004126506381281E-15, 3.09722235760629949E-12, 
                         4.60117603486559168E-10, 2.13941944795610622E-08, 2.46764213457981401E-07, 
                         2.73422068011878881E-06, 3.57293481989753322E-05, 0.000275242141167851312, 
                         0.000818953927502267349, 0.00231134524035220713,  0.00315544626918755127, 
                         0.015673473751851151,    0.0452736854651503914,   0.0923647267169863534, 
                         0.148070831155215854,    0.191760115888044341,    
                         0.000514894508069213769, 
                         0.191760115888044341,    0.148070831155215854,    
                         0.0923647267169863534,   0.0452736854651503914,   0.015673473751851151,    
                         0.00315544626918755127,  0.00231134524035220713,  0.000818953927502267349,
                         0.000275242141167851312, 3.57293481989753322E-05, 2.73422068011878881E-06,     
                         2.46764213457981401E-07, 2.13941944795610622E-08, 4.60117603486559168E-10,     
                         3.09722235760629949E-12, 5.45004126506381281E-15, 1.05413265823340136E-18 };



void getCC ( int n, int *nq, double **x, double **w ) {

  if ((n-1)%2 != 0) std::cout<<n<<std::endl;
  assert((n-1)%2 == 0);

  if ( *x != NULL ) free (*x);
  if ( *w != NULL ) free (*w);
  (*x) = (double *) malloc (n*sizeof(double));
  (*w) = (double *) malloc (n*sizeof(double));

  if ( n == 1 ) {
    (*x)[0] = 0.0; (*w)[0]=2.0; *nq = 1;
  }
  else if (( n == 3 ) || (n==5) || (n==9) || (n==17) || (n==33) || (n==65)){
    double dpi = 4.0*atan(1.0);
    int nm1 = n-1;
    for (int i=0; i<n; i++) { 
      (*x)[i] = cos(i*dpi/nm1);
      if ( fabs((*x)[i])<1.e-15) (*x)[i] =0.0;
    }

    (*w)[0  ] = 1.0/(nm1*nm1-1.0);
    (*w)[nm1] = (*w)[0];
    for (int i=1; i<nm1; i++) {
      (*w)[i] = 1.0+cos(i*dpi)/(1.0-nm1*nm1);
      for (int j = 1; j<=nm1/2-1; j++) 
	(*w)[i] += (2.0/(1.0-4.*j*j))*cos(2.*i*j*dpi/nm1);
      (*w)[i] *= 2.0/nm1;
    }

    *nq=n;
  }
  else
  {
    std::cout<<"Error in getCC(): n be one of {1,3,5,9,17,33,65}: "<<n<<std::endl;
    exit (1) ;
  }

  return;

}

void getGKPunif ( int n, int *nq, double **x, double **w ) {

  if ( n == 1 ) {
    *x = x1; *w=w1; *nq = 1;
  }
  else if ( n == 3 ) {
    *x = x3; *w=w3; *nq = 3;
  }
  else if ( n == 7 ) {
    *x = x7; *w=w7; *nq = 7;
  }
  else if ( n == 15 ) {
    *x = x15; *w=w15; *nq = 15;
  }
  else if ( n == 31 )
  {
    *x = x31; *w=w31; *nq = 31;
  }
  else if ( n == 63 )
  {
    *x = x63; *w=w63; *nq = 63;
  }
  else
  {
    std::cout<<"Error in getGKPunif(): n be one of {1,3,7,15,31,63}: "<<n<<std::endl;
    exit (1) ;
  }

  return;

}

void getGKPnorm ( int n, int *nq, double **x, double **w ) {
  if ( n == 1 ) {
    *x = xn1; *w=wn1; *nq = 1;
  }
  else if ( n == 3 ) {
    *x = xn3; *w=wn3; *nq = 3;
  }
  else if ( n == 9 ) {
    *x = xn9; *w=wn9; *nq = 9;
  }
  else if ( n == 19 ) {
    *x = xn19; *w=wn19; *nq = 19;
  }
  else if ( n == 35 )
  {
    *x = xn35; *w=wn35; *nq = 35;
  }
  else
  {
    std::cout<<"Error in getGKPnorm(): n be one of {1,3,9,19,35}: "<<n<<std::endl;
    exit (1) ;
  }

  return;

}

int getOrderCC ( int lev ) {

  if ( lev == 1 ) return (1);
  if ( lev <= 7 ) return ( (int) pow(2,lev-1) +1 );
  
  std::cout<<"Error in getOrderCC() : Level should be between 1...7: "<<lev<<endl;
  exit (-1) ;
  
}

int getOrderGKPunif ( int lev ) {

  if ( lev==1  ) return (1 );
  if ( lev<=3  ) return (3 );
  if ( lev<=6  ) return (7 );
  if ( lev<=12 ) return (15);
  if ( lev<=24 ) return (31);
  if ( lev<=48 ) return (63);
  
  std::cout<<"Error in getOrderGKPunif() : Level should be between 1...48: "<<lev<<endl;
  exit (-1) ;
  
}

int getOrderGKPnorm ( int lev ) {

  if ( lev==1  ) return (1 );
  if ( lev<=3  ) return (3 );
  if ( lev<=8  ) return (9 );
  if ( lev<=14 ) return (19);
  if ( lev<=25 ) return (35);
  
  std::cout<<"Error in getOrderGKPnorm() : Level should be between 1...48: "<<lev<<endl;
  exit (-1) ;
  
}

void getCompNintoDim(int n, int dim, int *nelem, int **plist) {
    
  *nelem = choose ( n + dim - 1, n );
  if (*plist != NULL) free(*plist);
  *plist = (int *) malloc((*nelem)*dim*sizeof(int));
  for (int i=0; i<(*nelem)*dim; i++) (*plist)[i]=0;

  (*plist)[0] = n;
  int t = n;
  int h = -1;
  int j = 0;

  while ( (*plist)[j*dim+dim-1]<n ) {
    j = j+1;
    for (int i = j*dim; i<(j+1)*dim; i++) (*plist)[i] = (*plist)[i-dim];
    if (t != 1) h = -1;
    h = h+1;
    t = (*plist)[j*dim+h];
    (*plist)[j*dim+h  ] = 0;
    (*plist)[j*dim    ] = t-1;
    (*plist)[j*dim+h+1] += 1;
  }

  for (int i=0; i<(*nelem)*dim; i++) (*plist)[i]+=1;


  return ;

}        

int getSpgSize ( int getOrder ( int ), int dim, int lev ) {

  /* 1D numbers */
  int *nq1D = (int *) malloc(lev*sizeof(int));
    for ( int i = 1; i<lev+1; i++ ) nq1D[i-1] = getOrder(i) ;

  /* Loop through all levels */
  int levMin = MAX ( 0, lev - dim );
  int spgSize = 0;

  for ( int ilev = levMin; ilev < lev; ilev++ ) {

    int *qlist=NULL, nelem;
    getCompNintoDim( ilev, dim, &nelem, &qlist );

    for ( int i=0; i<nelem;  i++ ) {
      int nqp = nq1D[qlist[i*dim]-1];
      for ( int j=1; j<dim; j++ ) nqp *= nq1D[qlist[i*dim+j]-1];
      spgSize += nqp;
    }

    if (qlist != NULL ) free(qlist);

  } 

  if (nq1D != NULL) free ( nq1D );

  return spgSize;

}

void getSpgQW ( void get1DQW ( int , int *, double **, double** ), int getOrder ( int ),
		int dim, int lev, int *nqpts, double **qpts, double **w ) {

#ifdef DEBUG
  std::cout<<"In the function:"<<dim<<" "<<lev<<std::endl;
#endif

  /* Initial estimate for number of quad points */
  (*nqpts) = getSpgSize ( getOrder, dim, lev );
#ifdef DEBUG
  std::cout<<(*nqpts)<<std::endl;
#endif

  if (*qpts != NULL) free(*qpts);
  *qpts = (double *) malloc((*nqpts)*dim*sizeof(double));
  if (*w != NULL) free(*w);
  *w = (double *) malloc((*nqpts)*sizeof(double));

  /* zero-out memory */
  for (int i=0; i<(*nqpts)*dim; i++ ) (*qpts)[i] = 0.0;
  for (int i=0; i<(*nqpts);     i++ ) (*w)[i]    = 1.0;

#ifdef DEBUG
  std::cout<<"Done with memory"<<std::endl;
#endif

  /* Create arrays to hold pointers to 1D quadrature rules */
  double **x1D = (double **) malloc(dim*sizeof(double *));
  double **w1D = (double **) malloc(dim*sizeof(double *));
  int    *n1D  = (int *)     malloc(dim*sizeof(int)     );

  for (int i=0; i<dim; i++) x1D[i] = NULL;
  for (int i=0; i<dim; i++) w1D[i] = NULL;

#ifdef DEBUG
  std::cout<<"Done with memory 2"<<std::endl;
#endif

  /* Loop through all levels */
  int levMin = MAX ( 0, lev - dim );
  int spgSize = 0;

  for ( int ilev = levMin; ilev < lev; ilev++ ) {

    int *qlist=NULL, nelem;
    getCompNintoDim( ilev, dim, &nelem, &qlist );

#ifdef DEBUG
    for (int ielem=0;ielem<nelem;ielem++) {
      for (int idim=0;idim<dim; idim++)
        std::cout<<qlist[ielem*dim+idim]<<" ";
      std::cout<<std::endl;
    }
#endif
    
    for ( int i=0; i<nelem;  i++ ) {

      /* Get pointers to 1D rules */
      for ( int j=0; j<dim; j++ ) 
        get1DQW ( getOrder(qlist[i*dim+j]), &n1D[j], &x1D[j], &w1D[j] );
      
      double qfac = choose ( dim - 1, dim + ilev - lev );
      if ( (lev-1-ilev)%2==1) qfac = -qfac;
      getTensorProd(dim,*qpts,*w,&spgSize,n1D,x1D,w1D,qfac);

    }

#ifdef DEBUG
    std::cout<<"Done with ilev:"<<ilev<<":"<<spgSize<<std::endl;
#endif

    if (qlist != NULL ) free(qlist);

  } 
  assert(spgSize == (*nqpts));
  sortSpg ( dim, spgSize, *qpts, *w );

  /* delete duplicates */
  int count = 0;
  for ( int i = 1; i < (*nqpts); i++ ) {
    bool isDup = true;
    for ( int j = 0; j < dim; j++ ) {
      if ( (*qpts)[count*dim+j] != (*qpts)[i*dim+j] ) {
        isDup = false;
        break;
      }
    }
    if ( isDup ) {
      (*w)[count] += (*w)[i];
    }
    else {
      count += 1;
      (*w)[count] = (*w)[i];
      for ( int j = 0; j < dim; j++ )
        (*qpts)[count*dim+j] = (*qpts)[i*dim+j];
    }
  }
  *nqpts = count + 1;

  return;

}

void getSpgAnisQW ( void get1DQW ( int , int *, double **, double** ), int getOrder ( int ),
		   int dim, int *levList, int *nqpts, double **qpts, double **w ) {

#ifdef DEBUG
  std::cout<<"In the function:"<<dim<<" "<<lev<<std::endl;
#endif

  int lev = 1;
  for ( int j=0; j<dim; j++ ) 
        if (levList[j]>lev) lev=levList[j];

  /* Initial estimate for number of quad points */
  (*nqpts) = getSpgSize ( getOrder, dim, lev );
#ifdef DEBUG
  std::cout<<(*nqpts)<<std::endl;
#endif

  if (*qpts != NULL) free(*qpts);
  *qpts = (double *) malloc((*nqpts)*dim*sizeof(double));
  if (*w != NULL) free(*w);
  *w = (double *) malloc((*nqpts)*sizeof(double));

  /* zero-out memory */
  for (int i=0; i<(*nqpts)*dim; i++ ) (*qpts)[i] = 0.0;
  for (int i=0; i<(*nqpts);     i++ ) (*w)[i]    = 1.0;

#ifdef DEBUG
  std::cout<<"Done with memory"<<std::endl;
#endif

  /* Create arrays to hold pointers to 1D quadrature rules */
  double **x1D = (double **) malloc(dim*sizeof(double *));
  double **w1D = (double **) malloc(dim*sizeof(double *));
  int    *n1D  = (int *)     malloc(dim*sizeof(int)     );

  for (int i=0; i<dim; i++) x1D[i] = NULL;
  for (int i=0; i<dim; i++) w1D[i] = NULL;

#ifdef DEBUG
  std::cout<<"Done with memory 2"<<std::endl;
#endif

  /* Loop through all levels */
  int levMin = MAX ( 0, lev - dim );
  int spgSize = 0;

  for ( int ilev = levMin; ilev < lev; ilev++ ) {

    int *qlist=NULL, nelem;
    getCompNintoDim( ilev, dim, &nelem, &qlist );

#ifdef DEBUG
    for (int ielem=0;ielem<nelem;ielem++) {
      for (int idim=0;idim<dim; idim++)
        std::cout<<qlist[ielem*dim+idim]<<" ";
      std::cout<<std::endl;
    }
#endif
    
    for ( int i=0; i<nelem;  i++ ) {

      bool goodElem=true;
      for ( int j=0; j<dim; j++ ) 
        if (qlist[i*dim+j]>levList[j]) {
	  goodElem = false;
	  std::cout<<i<<" "<<j<<": "<<dim<<":"<<qlist[i*dim+j]<<" "<<levList[j]<<std::endl;
	}

      if ( !goodElem ) continue;
      /* Get pointers to 1D rules */
      for ( int j=0; j<dim; j++ ) 
        get1DQW ( getOrder(qlist[i*dim+j]), &n1D[j], &x1D[j], &w1D[j] );
      
      double qfac = choose ( dim - 1, dim + ilev - lev );
      if ( (lev-1-ilev)%2==1) qfac = -qfac;
      getTensorProd(dim,*qpts,*w,&spgSize,n1D,x1D,w1D,qfac);

    }

#ifdef DEBUG
    std::cout<<"Done with ilev:"<<ilev<<":"<<spgSize<<std::endl;
#endif

    if (qlist != NULL ) free(qlist);

  } 
  //assert(spgSize == (*nqpts));
  std::cout<<spgSize<<std::endl;
  sortSpg ( dim, spgSize, *qpts, *w );

  /* delete duplicates */
  int count = 0;
  for ( int i = 1; i < (*nqpts); i++ ) {
    bool isDup = true;
    for ( int j = 0; j < dim; j++ ) {
      if ( (*qpts)[count*dim+j] != (*qpts)[i*dim+j] ) {
        isDup = false;
        break;
      }
    }
    if ( isDup ) {
      (*w)[count] += (*w)[i];
    }
    else {
      count += 1;
      (*w)[count] = (*w)[i];
      for ( int j = 0; j < dim; j++ )
        (*qpts)[count*dim+j] = (*qpts)[i*dim+j];
    }
  }
  *nqpts = count + 1;

  return;

}


void sortSpg ( int dim, int spgSize, double *qpts, double *w ) {

  assert(dim>0);
  assert(spgSize>1);

  int isgn=0, i1=0, j1=0, index=0;

  do {

    heap_ext_(&spgSize,&isgn,&i1,&j1,&index);
    if (index < 0) {
      isgn = 0;
      for ( int j = 0; j < dim; j++ ) {
        if ( qpts[(i1-1)*dim+j] < qpts[(j1-1)*dim+j] ) {
          isgn = -1; break;
        }
        else if ( qpts[(j1-1)*dim+j] < qpts[(i1-1)*dim+j] ) {
          isgn = 1;
          break;
        }
      }
    }

    if (index > 0) {
      double dtmp ;
      for ( int j = 0; j < dim; j++ ) {
        double dtmp = qpts[(i1-1)*dim+j];
        qpts[(i1-1)*dim+j] = qpts[(j1-1)*dim+j];
        qpts[(j1-1)*dim+j] = dtmp;
      }
      dtmp = w[i1-1];
      w[i1-1] = w[j1-1];
      w[j1-1] = dtmp;
    }

  } while (index != 0);


  return;

}

void getTensorProd(int dim, double *qpts, double *w, int *spgSize, int *n1D,
                   double **x1D, double **w1D, double qfac) {

  int n1 = 1, n2 = 1;
  for (int i=1; i<dim; i++ ) n2 *= n1D[i];

#ifdef DEBUG
  for (int idim=0;idim<dim;idim++) {
    std::cout<<"getTensorProd: dim="<<idim<<": "<<n1D[idim]<<std::endl;
    for (int i=0;i<n1D[idim]; i++)
        std::cout<<x1D[idim][i]<<" "<<w1D[idim][i]<<std::endl;
  }
#endif
  

  int ist = (*spgSize);
  for ( int k=0; k<dim; k++ ) {
    int idx = 0;
    for (int i=0; i<n1; i++ )
      for (int l=0; l<n1D[k]; l++ )
        for (int j=0; j<n2; j++) {
          int idxQ = ((*spgSize)+idx)*dim+k;
          qpts[idxQ] = x1D[k][l];
          w[(*spgSize)+idx] *= w1D[k][l];
          idx ++;
        }
    if (k<dim-1) {
      n1 *= n1D[k];
      n2 /= n1D[k+1];
    }
  }

  /* multiply by appropriate weights */
  n1 = n1D[0];
  for (int i=1; i<dim; i++ ) n1 *= n1D[i];
  for (int i=(*spgSize); i<(*spgSize)+n1; i++) w[i]*=qfac;
  (*spgSize) += n1;

  return ;

}
