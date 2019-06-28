#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <assert.h>
#include <algorithm>
#include <cmath>
#ifdef USE_HDF5
#include "H5Cpp.h"
using namespace H5;
#endif
typedef std::vector<double>  RealVector;
typedef std::vector<int>      IntVector;

char *getCmdOption(char **begin, char **end, const std::string &option) {

  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    return *itr;

  return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}

double tmcmc(RealVector &rngs, double gm, int nspl, RealVector dtvec, int iseed,
          const int nProcs, int ndim, bool spc);

double tmcmc(RealVector &rngs, double gm, int nspl, RealVector dtvec, int iseed,
          const int nProcs, int ndim, bool spc, double cv, double a, double b,
          bool mala, double tauScale, double betaThres, int MFactor, bool basis, int CATSteps);

int main(int argc, char *argv[]) {

  // Default Arguments
  int iseed  = 2014;
  int ndim   = 2;
  int nprocs = 2;
  int nspl   = 50000;
  double gamma = 0.2;
  int maxniter = 1;
  double initDelta = 2;
  double deltaRate = 0.8;
  int spc = 0; // Clustering Analysis, boolean
  // Stress Test Purposes
  double cv = 1;
  double a = 1.0/9.0;
  double b = 8.0/9.0;
  bool mala = false;
  double tauScale = 0.05;
  double betaThres = 0.5;
  int MFactor = 1;
  bool basis = false;
  int CATSteps = 1;

  // command-line options
  if( cmdOptionExists(argv, argv+argc,"-s"))
    iseed=atoi(getCmdOption(argv,argv + argc,"-s"));
  if( cmdOptionExists(argv, argv+argc,"-d"))
    ndim=atoi(getCmdOption(argv,argv + argc,"-d"));
  if( cmdOptionExists(argv, argv+argc,"-p"))
    nprocs=atoi(getCmdOption(argv,argv + argc,"-p"));
  if( cmdOptionExists(argv, argv+argc,"-n"))
    nspl=atoi(getCmdOption(argv,argv + argc,"-n"));
  if( cmdOptionExists(argv, argv+argc,"-g"))
    gamma=atof(getCmdOption(argv,argv + argc,"-g"));
  if( cmdOptionExists(argv, argv+argc,"-i"))
    maxniter=atoi(getCmdOption(argv, argv + argc, "-i"));
  if( cmdOptionExists(argv, argv+argc,"-t"))
    initDelta=atof(getCmdOption(argv, argv + argc, "-t"));
  if( cmdOptionExists(argv, argv+argc,"-r"))
    deltaRate=atof(getCmdOption(argv, argv + argc, "-r"));
  if( cmdOptionExists(argv, argv+argc, "-c"))
    spc=atoi(getCmdOption(argv, argv + argc, "-c"));
  if( cmdOptionExists(argv, argv+argc, "-v"))
    cv=atof(getCmdOption(argv, argv + argc, "-v"));
  if( cmdOptionExists(argv, argv+argc, "-a"))
    a=(double)1.0/(double)9.0 * atof(getCmdOption(argv, argv + argc, "-a"));
  if( cmdOptionExists(argv, argv+argc, "-b"))
    b=(double)8.0/(double)9.0 * atof(getCmdOption(argv, argv + argc, "-b"));
  if( cmdOptionExists(argv, argv+argc,"-m"))
    mala=atoi(getCmdOption(argv,argv + argc,"-m"));
  if( cmdOptionExists(argv, argv+argc,"-B"))
    betaThres=atof(getCmdOption(argv,argv + argc,"-B"));
  if( cmdOptionExists(argv, argv+argc,"-M"))
    MFactor=atof(getCmdOption(argv,argv + argc,"-M"));
  if( cmdOptionExists(argv, argv+argc,"-l"))
    basis=atoi(getCmdOption(argv,argv + argc,"-l"));
  if( cmdOptionExists(argv, argv+argc,"-C"))
    CATSteps=atoi(getCmdOption(argv,argv + argc,"-C"));

  std::cout<<"No. of samples   : "<<nspl <<std::endl;
  std::cout<<"No. of processors: "<<nprocs<<std::endl;
  std::cout<<"Gamma factor     : "<<gamma <<std::endl;
  std::cout<<"Beta Threshold   : "<<betaThres<<std::endl;

  RealVector rngs;
  std::ifstream myfile;
  double dtmp;
  myfile.open("rngs.dat", std::ios_base::in);
  while (myfile >> dtmp) rngs.push_back(dtmp);
  myfile.close();

  std::cout<<"Parameter ranges"<<std::endl;
  for (int i=0; i<(rngs.size()/2); i++) {
    std::cout<<"Par "<<i+1<<": "<<rngs[2*i]<<"<->"<<rngs[2*i+1]<<std::endl;
  }



  #ifdef USE_HDF5
  // HDF5 Attribute Initialization
  H5File h5file("RunResults.h5", H5F_ACC_TRUNC);
  hsize_t attrDims[1] = {1};
  hsize_t rngsDims[2] = {ndim, 2};
  DataSpace attrDataSpace = DataSpace (1, attrDims);
  DataSpace rngsDataSpace = DataSpace (2, rngsDims);

  Attribute initDeltaAttr = h5file.createAttribute("InitDelta",
              PredType::NATIVE_DOUBLE, attrDataSpace);
  Attribute deltaRateAttr = h5file.createAttribute("DeltaRate",
              PredType::NATIVE_INT, attrDataSpace);
  Attribute nsplAttr = h5file.createAttribute("Number of Samples",
              PredType::NATIVE_INT, attrDataSpace);
  Attribute ndimAttr = h5file.createAttribute("Number of Dimensions",
              PredType::NATIVE_INT, attrDataSpace);
  Attribute seedAttr = h5file.createAttribute("Random Seed",
              PredType::NATIVE_INT, attrDataSpace);
  Attribute nstagesAttr = h5file.createAttribute("Number of Stages",
              PredType::NATIVE_INT, attrDataSpace);
  Attribute spcAttr = h5file.createAttribute("Spectral Clustering",
              PredType::NATIVE_INT, attrDataSpace);
  DataSet rngsDataset = h5file.createDataSet("Ranges",
              PredType::NATIVE_DOUBLE, rngsDataSpace);

  int iseedAttrData[1]  = {iseed};
  int ndimAttrData[1]   = {ndim};
  int nsplAttrData[1]   = {nspl};
  int maxniterAttrData[1] = {maxniter};
  double initDeltaAttrData[1] = {initDelta};
  double deltaRateAttrData[1] = {deltaRate};
  int spcAttrData[1] = {spc};
  double rngsData[ndim][2];
  for (int i=0; i<ndim; i++) {
    rngsData[i][0] = rngs[2*i];
    rngsData[i][1] = rngs[2*i + 1];
  }

  initDeltaAttr.write(PredType::NATIVE_DOUBLE, initDeltaAttrData);
  deltaRateAttr.write(PredType::NATIVE_DOUBLE, deltaRateAttrData);
  nsplAttr.write(PredType::NATIVE_INT, nsplAttrData);
  ndimAttr.write(PredType::NATIVE_INT, ndimAttrData);
  seedAttr.write(PredType::NATIVE_INT, iseedAttrData);
  nstagesAttr.write(PredType::NATIVE_INT, maxniterAttrData);
  spcAttr.write(PredType::NATIVE_INT, spcAttrData);
  rngsDataset.write(rngsData, PredType::NATIVE_DOUBLE);
  rngsDataset.close();
  #endif
  if (maxniter > 1) {
    std::ofstream deltaRateFile("deltaRate.dat");
    deltaRateFile.clear();
    deltaRateFile << deltaRate;
    deltaRateFile.close();
  }
  RealVector dts; // Obsolete, Dbeta changes are adaptive now.
  std::ofstream evidFile("Evidence.dat");

  // Run TMCMC (Repeat for Manifolds)
  for (size_t i = 0; i < maxniter; ++i) {


    #ifdef USE_HDF5
    // Track stage with HDF5
    std::string stageName = "Stage" + std::to_string(i);
    Group stage = h5file.createGroup(stageName.c_str());
    // Track Delta
    Attribute deltaAttr = stage.createAttribute("Current Delta",
                        PredType::NATIVE_DOUBLE, attrDataSpace);
    double deltaAttrData[1] = {initDelta * pow(deltaRate,i)};
    deltaAttr.write(PredType::NATIVE_DOUBLE, deltaAttrData);
    // stage.close();
    // h5file.close();
    #endif

    // Track Stage
    std::ofstream TMCMCiter("TMCMCiter.dat");
    TMCMCiter.clear();
    TMCMCiter << i;
    TMCMCiter.close();
    // Track Current Delta
    std::ofstream deltaFile("delta.dat");
    deltaFile.clear();
    deltaFile << std::setprecision(18) << initDelta * pow(deltaRate, i);
    deltaFile.close();

    // Run TMCMC, get evidence
    double evid = tmcmc(rngs, gamma, nspl, dts, iseed, nprocs, ndim, spc,
                        cv, a, b, mala, tauScale, betaThres, MFactor, basis, CATSteps);

    #ifdef USE_HDF5
    stage = h5file.openGroup(stageName.c_str());
    Attribute evidAttr = stage.createAttribute("Evidence",
                          PredType::NATIVE_DOUBLE, attrDataSpace);
    double evidAttrData[1] = {evid};
    evidAttr.write(PredType::NATIVE_DOUBLE, evidAttrData);
    stage.close();
    #else
    evidFile << std::setprecision(18) << evid << std::endl;
    // Move residual files out of the way if necessary.
    std::string filename;
    for (size_t f = 1000; f > 0; --f) {
      filename = "samples.dat." + std::to_string(f);

      std::ifstream sampleFile(filename);

      if (sampleFile.is_open()) {
        break;
      }
    }

    std::ifstream moveFile("move.sh");
    if (moveFile.is_open()) {
      std::string moveStr = "./move.sh Stage" + std::to_string(i) +
      " " + filename;
      system(moveStr.c_str());
    }
    #endif
  }

  #ifdef USE_HDF5
  std::string rmstr = "./hdf5clean.sh";
  std::cout << "CLEANUP";
  system(rmstr.c_str());
  #else
  evidFile.close();
  #endif

  return (0);
}
