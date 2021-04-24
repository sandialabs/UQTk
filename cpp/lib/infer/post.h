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
/// \file post.h
/// \author K. Sargsyan  2016 -
/// \brief Header for the posterior computation class

#ifndef POST_H_SEEN
#define POST_H_SEEN

#include "Array1D.h"
#include "Array2D.h"
#include "mrv.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <sstream>

using namespace std; // needed for python string conversion

/// \class	Post: class for posterior evaluation with various likelihood and prior options

class Post {
public:

  /// \brief Constructor
  Post();
  /// \brief Destructor
  ~Post() {}

  /// \brief Set the x- and y-data
  void setData(Array2D<double>& xdata,Array2D<double>& ydata);
  /// \brief Set the x- and y-data, for the case when each datapoint have different number of measurements
  void setData(Array2D<double>& xdata,Array1D<Array1D<double> >& ydata);
  /// \brief Set the magnitude of data noise
  void setDataNoise(Array1D<double>& sigma);
  /// \brief Indicate inference of data noise stdev
  void inferDataNoise();
  /// \brief Indicate inference of log of data noise stdev
  void inferLogDataNoise();
  /// \brief Get data noise, whether inferred or fixed
  Array1D<double> dataSigma(double m_last);
  /// \brief Set a pointer to the forward model f(p,x)
  void setModel(Array1D< Array2D<double> (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *) > forwardFuncs, Array2D<double>& fixindnom, void* funcInfo);
  /// \brief Set model input parameters' randomization scheme
  void setModelRVinput(int pdim, int order, Array1D<int>& rndInd,string pdfType,string pcType);
  /// \brief Get the dimensionailty of the posterior function
  int getChainDim();
  /// \brief Set the prior type and its parameters
  void setPrior(string priorType, double priora, double priorb);
  /// \brief Evaluate log-prior
  double evalLogPrior(Array1D<double>& m);
  /// \brief Extract parameter PC coefficients from a posterior input
  Array2D<double> getParamPCcf(Array1D<double>& m);
  /// \brief Sample model parameters given posterior input
  Array2D<double> samParam(Array1D<double>& m, int ns);
  /// \brief Get moments of parameters given posterior input
  void momParam(Array1D<double>& m, Array1D<double>& parMean, Array1D<double>& parVar, bool covFlag, Array2D<double>& parCov);
  /// \brief Sample forward function at a given grid for given posterior input
  Array2D<double> samForwardFcn(Array2D<double> (*forwardFunc)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*),Array1D<double>& m, Array2D<double>& xgrid, int ns);
  /// \brief Get moments of forward function at a given grid for given posterior input
  void momForwardFcn(Array2D<double> (*forwardFunc)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*),Array1D<double>& m, Array2D<double>& xgrid, Array1D<double>& fcnMean, Array1D<double>& fcnVar, bool covflag, Array2D<double>& fcnCov);
  /// \brief Get moments of composite forward function at a given grid for given posterior input
  void momForwardFcn(Array1D<double>& m, Array2D<double>& xgrid, Array1D<double>& fcnMean, Array1D<double>& fcnVar, bool covflag, Array2D<double>& fcnCov);
  /// \brief Dummy evaluation of log-likelihood
  virtual double evalLogLik(Array1D<double>& m){return 0;};


protected:
  /// \brief xdata
  Array2D<double> xData_;
  /// \brief ydata
  Array1D<Array1D<double>> yData_;
  /// \brief ydata averaged per measurement
  Array1D<double> yDatam_;
  /// \brief Number of data points
  int nData_;
  /// \brief Number of samples at each input
  Array1D<int> nEachs_;
  /// \brief Dimensionality of x-space
  int xDim_;
  /// \brief Dimensionality of parameter space (p-space)
  int pDim_;
  /// \brief Dimensionality of posterior input
  int chDim_;
  /// \brief Flag for data noise inference
  bool inferDataNoise_;
  /// \brief Flag to check if data noise logarithm is used
  bool dataNoiseLogFlag_;
  /// \brief Data noise stdev
  Array1D<double> dataNoiseSig_;
  /// \brief Pointer to the forward function f(p,x)
  //Array2D<double> (*forwardFcn_)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void*);
  Array1D< Array2D<double> (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *) > forwardFcns_;
  /// \brief Auxiliary information for function evaluation
  void* funcinfo_;
  /// \brief Number of extra inferred parameters, such as data noise or Koh variance
  int extraInferredParams_;
  /// \brief Number of categories
  int ncat_;

  /// \brief Pointer to a multivariate PC RV object
  Mrv * Mrv_;
  /// \brief Indices of randomized inputs
  Array1D<int> rndInd_;
  /// \brief Indices and nominal values for fixed inputs
  Array2D<double> fixIndNom_;
  /// \brief Lower and upper bounds on parameters
  Array1D<double> lower_,upper_;
  /// \brief Input parameter PDF type
  string pdfType_;
  /// \brief PC type parameter for the r.v.
  string rvpcType_;
  /// \brief Prior type
  string priorType_;
  /// \brief Prior parameter #1
  double priora_;
  /// \brief Prior parameter #2
  double priorb_;

 private:

  /// \brief Verbosity level
  int verbosity_;

};

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

/// \class Lik_Full
/// \brief Derived class for full likelihood
class Lik_Full: public Post {
public:
  /// \brief Constructor given KDE bandwidth and sample size
  Lik_Full(double bdw,int nsam){this->bdw_=bdw; this->nsam_=nsam; return;}
  /// \brief Destructor
  ~Lik_Full(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);

  private:
    /// \brief KDE bandwidth
    double bdw_;
    /// \brief KDE sample size
    int nsam_;
};

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

/// \class Lik_Marg
/// \brief Derived class for marginal likelihood
class Lik_Marg: public Post {
public:
  /// \brief Constructor given KDE bandwidth and sample size
  Lik_Marg(double bdw,int nsam){this->bdw_=bdw; this->nsam_=nsam; return;}
  /// \brief Destructor
  ~Lik_Marg(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);

  private:
    /// \brief KDE bandwidth
    double bdw_;
    /// \brief KDE sample size
    int nsam_;

};

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

/// \class Lik_MVN
/// \brief Derived class for mvn likelihood
class Lik_MVN: public Post {
public:
  /// \brief Constructor given fiagonal nugget
  Lik_MVN(double nugget){this->nugget_=nugget; return;}
  /// \brief Destructor
  ~Lik_MVN(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);

  private:
    /// \brief Nugget size
    double nugget_;

};

/*******************************************************************/

/// \class Lik_GausMarg
/// \brief Derived class for gaussian-marginal likelihood
class Lik_GausMarg: public Post {
public:
  /// \brief Constructor
  Lik_GausMarg(){return;}
  /// \brief Destructor
  ~Lik_GausMarg(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);

};

/*******************************************************************/

/// \class Lik_GausMargD
/// \brief Derived class for gaussian-marginal likelihood with discrete parameter
class Lik_GausMargD: public Post {
public:
  /// \brief Constructor
  Lik_GausMargD(){return;}
  /// \brief Destructor
  ~Lik_GausMargD(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);

};

/*******************************************************************/

/// \class Lik_ABC
/// \brief Derived class for ABC likelihood
class Lik_ABC: public Post {
public:
  /// \brief Constructor given ABC epsilon
  Lik_ABC(double eps){this->abceps_=eps; return;}
  /// \brief Destructor
  ~Lik_ABC(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);

private:
  /// \brief ABC epsilon
  double abceps_;

};

/*******************************************************************/

/// \class Lik_ABCm
/// \brief Derived class for ABC-mean likelihood
class Lik_ABCm: public Post {
public:
  /// \brief Constructor given ABC epsilon
  Lik_ABCm(double eps){this->abceps_=eps; return;}
    /// \brief Destructor
  ~Lik_ABCm(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);

private:
  /// \brief ABC epsilon
  double abceps_;

};

/*******************************************************************/

/// \class Lik_Koh
/// \brief Derived class for Kennedy-O'Hagan likelihood
class Lik_Koh: public Post {
public:
  /// \brief Constructor given correlation length
  Lik_Koh(double corLength){  this->extraInferredParams_=1;this->corLength_=corLength; return;}
  /// \brief Destructor
  ~Lik_Koh(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);
private:
  double corLength_;

};

/*******************************************************************/

/// \class Lik_Classical
/// \brief Derived class for classical likelihood
class Lik_Classical: public Post {
public:
  /// \brief Constructor
  Lik_Classical(){ return;}
  /// \brief Destructor
  ~Lik_Classical(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);


};

/*******************************************************************/
/// \class Lik_Eov
/// \brief Derived class for error-in-variable likelihood
class Lik_Eov: public Post {
public:
  /// \brief Constructor
  Lik_Eov(){ return;}
  /// \brief Destructor
  ~Lik_Eov(){};

  /// \brief Evaluate log-likelihood
  double evalLogLik(Array1D<double>& m);


};

/*******************************************************************/

#endif /* POST_H_SEEN */
