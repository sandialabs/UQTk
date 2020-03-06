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
#ifndef USER_FUNCS_H_
#define USER_FUNCS_H_

// specification of the nominal values for parameters
void userSpecifyNominalParams(dataPosteriorInformation & dataPostInfo);
// specification of constraint information on parameters (e.g. statistics etc.), and adherence coefficients (ABC deltas)
void userDefineConstraints(dataPosteriorInformation & dataPostInfo);
// specification of (noisy) data, e.g. read from file
void userDefineData(dataPosteriorInformation &dataPostInfo);
//  specification of the fitting model
void userRunModel(Array1D<double> &modelDatay, Array1D<double> & modelDatax, Array1D<double> &parameters, Array1D<double> &hyperparameters);
// specification of user-defined parameter posterior (likelihood x prior). Runs the model and computes the likelihood
double userComputeParamLogPosterior(parameterPosteriorInformation * paramPostInfo, Array1D<double> parameters);
// specification of user-defined parameter likelihood function
double userComputeParamLogLikelihood(parameterPosteriorInformation * paramPostInfo, Array1D<double> modelDataOut, Array1D<double> parameters, Array1D<double> hyperparameters);
// specification of user defined statistics of the parameter posterior (e.g. moments)
void userComputeStatistics(Array1D<double> &parameterStatistics, Array1D<MCMC::chainstate> & parameterChain);

#endif  //USER_FUNCS_H_
