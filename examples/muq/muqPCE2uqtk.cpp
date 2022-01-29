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
//Standard includes
#include <Eigen/Core>
#include <Eigen/QR>
#include <iostream>

// MUQ includes
#include <MUQ/Modelling/ModPieceTemplates.h>
#include <MUQ/Modelling/CachedModPiece.h>
#include <MUQ/Approximation/pce/PolynomialChaosExpansion.h>
#include <MUQ/Approximation/smolyak/SmolyakPCEFactory.h>
#include <MUQ/Approximation/utilities/PolychaosXMLHelpers.h>
#include <MUQ/Utilities/LogConfig.h>
#include "dsfmt_add.h"

// UQTk and dependencies include
#define HAVEUQTk
#ifdef HAVEUQTk
  #include "XMLUtils.h"
  #include "tools.h"
  #include "mcmc.h"
  #include "arrayio.h"
  #include "arraytools.h"
  #include "PCSet.h"
#endif

// MUQ and standard namespaces
using namespace muq::Modelling;
using namespace muq::Approximation;
using namespace std;

// DEFINES
#ifndef NINP
  #define NINP   5
#endif
#ifndef PCETOL
  #define PCETOL 1.e-2
#endif
#ifndef RNSEED
  #define RNSEED 1984
#endif
#ifndef NL2EVAL
  #define NL2EVAL 10000
#endif

// Define the model
class ExampleModPiece : public OneInputNoDerivModPiece {
public:

  ExampleModPiece() : OneInputNoDerivModPiece(NINP, 1) {}
  virtual ~ExampleModPiece() = default;

private:

  /*
   * The implementation of the model
   */
  virtual Eigen::VectorXd EvaluateImpl(Eigen::VectorXd const& input) override
  {

    Eigen::VectorXd rval(1);

    ostringstream s;
    double meval;
    s << "./genzGP.py  " ;
    for (int i=0; i<NINP; i++)
      s << input(i) << "  " ;
    system(s.str().c_str());

    fstream myfile;
    myfile.open ("rval.txt");
    myfile >> meval;
    myfile.close();

    rval << meval ;

    // double dsum = 0.0;
    // double ai[]={0.001,0.01,0.01,0.1,1.0};
    // for (int i=0; i<NINP; i++)
    //   dsum += ai[i]*input(i)*input(i);
    // rval << exp(-dsum) ;

    return rval;

  }

};


int main(int argc, char **argv)
{

  //Initialization boilerplate for logging. Works regardless of whether GLOG is installed.
  muq::Utilities::InitializeLogging(argc, argv);

  // Create model
  auto trueModel = make_shared<ExampleModPiece>();

  // Create input xi's
  boost::property_tree::ptree ptree;
  for (int i=0; i<NINP; i++) {
    boost::property_tree::ptree subtree_legendre;
    subtree_legendre.put("index", i+1);
    subtree_legendre.put("type", "Legendre");
    subtree_legendre.put("count", 1);
    ptree.add_child("Variable", subtree_legendre);
  }
  //cout<<"Done creating tree of xi's"<<endl;

  //VariableCollection: automatically constructs the required quadrature rules and polynomials
  auto variableCollection = muq::Approximation::ConstructVariableCollection(ptree);

  // Cache model evaluations
  auto cachedTrueModel = make_shared<CachedModPiece>(trueModel);

  // Smolyak object: perform the pseudospectral construction and make a PCE for us.
  auto pceFactory = make_shared<SmolyakPCEFactory>(variableCollection, cachedTrueModel);

  // Make PCE: up to desired tolerance; there are other Start... methods
  std::shared_ptr<PolynomialChaosExpansion> pce = pceFactory->StartAdaptiveToTolerance(2, PCETOL);

  // Ane can re-adapt to another tolerance
  //pce = pceFactory->AdaptToTolerance(0.01);

  // Evaluate L2 error
/*  dsfmt_gv_init_gen_rand(RNSEED);
  Eigen::VectorXd testPoint(NINP);
  double l2sumNom=0.0, l2sumDen=0.0;
  double modval, pceval;
  for (int iev = 0; iev < NL2EVAL; iev++) {
    for (int i=0; i<NINP; i++)
      testPoint(i) = dsfmt_gv_genrand_urv_sm(-1.0,1.0);
    modval = trueModel->Evaluate(testPoint)(0);
    pceval = pce->Evaluate(testPoint)(0);
    l2sumNom += (modval-pceval)*(modval-pceval);
    l2sumDen += modval*modval;
  }
  cout << "Adaptive Tolerance/No. of evals/L2 error " << PCETOL << ","
       << cachedTrueModel->GetNumOfEvals() <<","
       << sqrt(l2sumNom/l2sumDen) << endl;
*/
  // multi-indices and PCE coefficients
  Eigen::MatrixXu midxMUQ = pce->GetMultiIndices();
  Eigen::MatrixXd coeffs  = pce->GetCoefficients();
  // print the size of the multiindex and coefficients

#ifdef VERBOSE
  cout << "mindex size: " << midxMUQ.rows() << " x " << midxMUQ.cols() << endl;
  cout << "coeffs size: " << coeffs.rows()  << " x " << coeffs.cols() << endl;
  for ( int i=0; i<(int) midxMUQ.rows(); i++) {
    for ( int j=0; j<(int) midxMUQ.cols(); j++)
      cout<<midxMUQ(i,j)<<" ";
    cout<<": "<<coeffs.transpose()(i)<<endl;
  }
#endif

#ifdef HAVEUQTk
  // Copy to Array classes
  Array2D<int> midxUQTk((int) midxMUQ.rows(),(int) midxMUQ.cols());
  Array1D<double> cUQTk((int) midxMUQ.rows());
  for (int i = 0; i < (int) midxMUQ.rows(); i++) {
    cUQTk(i) = coeffs(0,i);
    for (int j = 0; j < (int) midxMUQ.cols(); j++) {
      midxUQTk(i,j) = midxMUQ(i,j);
    }
  }

  // // Create PCSet class with custom m-index
  // PCSet myPCSet("NISPnoq",midxUQTk,"LU",0.0,1.0);

  // Array1D<double> rPoint(NINP);
  // double mxerr = 0.0;
  // for (int iev = 0; iev < 100; iev++) {
  //   for (int i=0; i<NINP; i++)
  //     testPoint(i) = rPoint(i) = dsfmt_gv_genrand_urv_sm(-1.0,1.0);
  //   double err = pce->Evaluate(testPoint)(0) - myPCSet.EvalPC(cUQTk,rPoint);
  //   if ( fabs(err) > mxerr ) mxerr = fabs(err);
  // }
  // cout << "Max abs diff between UQTk and MUQ: "<<mxerr<<endl;
#endif

}
