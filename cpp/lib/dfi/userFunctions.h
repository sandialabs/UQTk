#ifndef USER_FUNCS_H_
#define USER_FUNCS_H_

void userSpecifyNominalParams(dataPosteriorInformation & dataPostInfo);
void userDefineConstraints(dataPosteriorInformation & dataPostInfo);
void userDefineData(dataPosteriorInformation &dataPostInfo);
void userRunModel(Array1D<double> &modelDatay, Array1D<double> & modelDatax, Array1D<double> &parameters, Array1D<double> &hyperparameters);
double userComputeParamLogPosterior(parameterPosteriorInformation * paramPostInfo, Array1D<double> parameters);
double userComputeParamLogLikelihood(parameterPosteriorInformation * paramPostInfo, Array1D<double> modelDataOut, Array1D<double> parameters, Array1D<double> hyperparameters);
void userComputeStatistics(Array1D<double> &parameterStatistics, Array1D<MCMC::chainstate> & parameterChain);

#endif  //USER_FUNCS_H_
