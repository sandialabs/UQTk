/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.0
                          Copyright (2020) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

struct postAux
     {
       Array2D<double> data;
       Array1D<double> modelparams;
       Array1D<string> modelparamnames;
       Array1D<double> modelauxparams;
       Array1D<double> postparams;
       string noisetype;
       Array1D<string> priortype;
       Array1D<double> priorparam1;
       Array1D<double> priorparam2;
       Array1D<int> chainParamInd;
 }; 



/// \brief Evaluate the log of the posterior for a given set of model
/// and nuisance parameters
double LogPosterior(Array1D<double>& m, void *info);



double LogPosterior(Array1D<double>& m, void *info)
{
  double pi=4.*atan(1.);

  postAux* pinfo = (postAux*) info;

  double logprior=0;
  
  for(int ic=0;ic<m.Length();ic++){

    if(!strcmp(pinfo->priortype(ic).c_str(),"uniform")){
      double a=pinfo->priorparam1(ic);
      double b=pinfo->priorparam2(ic);
      if (a>=b)
        throw Tantrum("Prior bounds are incorrect!");
      if(m(ic)<b && m(ic)>a)
        logprior -= log(b-a);
      else
         return -1.e180;
    }

    else if(!strcmp(pinfo->priortype(ic).c_str(),"normal")){
      double mu=pinfo->priorparam1(ic);
      double sig=pinfo->priorparam2(ic);
    
      logprior -= 0.5*log(2.*pi*sig*sig);
      logprior -= 0.5*pow((m(ic)-mu)/sig,2.0);
    }
    else
      throw Tantrum("Only uniform or normal priors are implemented!");
  }


  Array2D<double> data=pinfo->data;
  Array1D<double> modelparams=pinfo->modelparams;
  // Set the parameter values for the forward run
  int chaindim=m.XSize();
  for(int ic=0;ic<chaindim;ic++)
     modelparams(pinfo->chainParamInd(ic))=m(ic);


  // Posterior parameter 
  // Either Proportionality constant between signal and noise for likelihood construction
  // or standard deviation itself
  double stdpar=pinfo->postparams(0);


  // Compare the model data with measurement data according to the noise model
  // to get the likelihood.
  double std;

  double sum=logprior;
  int nTot=data.XSize();
  for(int it=0; it < nTot; it++){

    double modelval=model(data(it,0),modelparams);


    if(!strcmp(pinfo->noisetype.c_str(),"const_stn")){
      std=stdpar*fabs(modelval);
    }
    else    if(!strcmp(pinfo->noisetype.c_str(),"const_stdev")){
      std=stdpar;
    }
    else
      throw Tantrum("Noise type is not recognized!");
 
    double var=std*std;

    double err=data(it,1)-modelval;
    sum  = sum - 0.5 * log(2*pi) - 0.5 * log(var)  - pow(err,2)/(2.0*var);
 
  }

  return sum;
 
}
