/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.3
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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file model_inf.cpp
/// \author K. Sargsyan 2015 -
/// \brief Command-line utility for model parameter inference

#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <map>
#include <iostream>
#include <string>
#include <math.h>

#include <getopt.h>


#include "func.h"
#include "post.h"
#include "mrv.h"
#include "inference.h"

#include "mcmc.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"

using namespace std;



/// default model type
#define MODELTYPE "linear" //"exp", "exp_quad", "prop", "prop_quad", "pcx","pc", "pcs", linear", "bb", "heat_transfer1", "heat_transfer2", "exp_sketch"
/// default likelihood type
#define LIKTYPE "classical" //"classical", "abc", "gausmarg", "marg", "mvn", "full", "koh"
/// default likelihood parameter of type double
#define LIKPARAM_DBL 0.01 // cov. nugget for MVN, cor. length for koh likelihood, KDE bandwidth for full and marg(but then Optimal is chosen!!), and epsilon for ABC likilohoods
/// default likelihood parameter of type int
#define LIKPARAM_INT 1000   // number of KDE samples for full and marg
/// default prior type
#define PRIORTYPE "uniform" //"uniform", "normal", "inverse", "wishart"
/// default prior parameter #1
#define PRIORA -DBL_MAX
/// default prior parameter #2
#define PRIORB DBL_MAX
/// default xfile
#define XFILE "xdata.dat"
/// default yfile
#define YFILE "ydata.dat"
/// default model parameter dimensionailty
#define PDIM 2
/// default parameter pdf order
#define ORDER 1
/// default parameter pdf type
#define PDFTYPE "pct"
/// default parameter pc for parameter pdf
#define PCTYPE "HG"
/// default datanoise
#define DATANOISE 0.1
/// default number of MCMC steps
#define NMCMC 10000 //if 0, only does optimization
/// default MCMC gamma (i.e. proposal size factor) for adaptive MCMC
#define MCMCGAMMA 0.1
/// default burn-in ratio
#define FBURN 10
/// default chain thinning
#define NSTEP 5

/// Displays information about this program
int usage(){
    // Note that all letters of alphabet are used, except -v (todo: omplement -v for verbosity levels)
    printf("This program infers model parameters given data.\n");
    printf("usage: model_inf [-h]  [-f<modeltype>] [-l<liktype>] [-w<likParam_dbl>] [-k<likParam_int>]\n");
    printf("                 [-i<priortype>] [-a<priora>] [-b<priorb>] [-j<chaininitfile] [-z]\n");
    printf("                 [-x<xfile>] [-y<yfile>] [-t<xgridfile>] \n");
    printf("                 [-e<datanoise>] \n");
    printf("                 [-d<pdim>] [-o<order>] [-r<rndindfile>] [-v<fixindnomfile>] [-s<pdftype>] [-c<pctype>] \n");
    printf("                 [-m<nmcmc>] [-g<mcmcgamma>] \n");
    printf("                 [-q<pgridfile>] [-p<pchainfile>] [-u<fburn>] [-n<nstep>]\n");
    printf(" -h                 : Print out this help message \n");
    printf(" -f <modeltype>     : Model type (default=%s) \n",MODELTYPE);
    printf(" -l <liktype>       : Likelihood type (default=%s) \n",LIKTYPE);
    printf(" -w <likParam_dbl>  : Likelihood parameter of type double (default=%lg) \n",LIKPARAM_DBL);
    printf(" -k <likParam_int>  : Likelihood parameter of type int (default=%d) \n",LIKPARAM_INT);
    printf(" -i <priortype>     : Prior type (default=%s) \n",PRIORTYPE);
    printf(" -a <priora>        : Prior parameter #1 (default=%lg) \n", PRIORA);
    printf(" -b <priorb>        : Prior parameter #2 (default=%lg) \n", PRIORB);
    printf(" -j <chaininitfile> : Chain initial state \n");
    printf(" -z                 : Whether to prepend optimization to MCMC \n");
    printf(" -x <xfile>         : Input x-data (default=%s) \n",XFILE);
    printf(" -y <yfile>         : Input y-data (default=%s) \n",YFILE);
    printf(" -t <xgridfile>     : x-grid where the model is evaluated \n");
    printf(" -d <pdim>          : Model parameter dimensionality (default=%d) \n",PDIM);
    printf(" -o <order>         : NISP order (default=%d) \n",ORDER);
    printf(" -r <rndindfile>    : Indices of randomized parameters \n");
    printf(" -v <fixindnomfile> : Indices and nominals of fixed parameters \n");
    printf(" -s <pdftype>       : Parameter pdf type (default=%s) \n",PDFTYPE);
    printf(" -c <pctype>        : PC type for parameter PDF (default=%s) \n",PCTYPE);
    printf(" -e <datanoise>     : Data noise parameter\n");
    printf(" -m <nmcmc>         : Number of MCMC samples (default=%d) \n",NMCMC);
    printf(" -g <mcmcgamma>     : MCMC parameter gamma (default=%lg) \n",MCMCGAMMA);
    printf(" -q <pgridfile>     : Parameter grid where the posterior is computed \n");
    printf(" -p <pchainfile>    : Parameters file for which the model is evaluated \n");
    printf(" -u <fburn>         : Ratio of burn-in MCMC samples (default=%d) \n",FBURN);
    printf(" -n <nstep>         : Thinning of MCMC samples (default=%d) \n",NSTEP);
    printf("=========================================================================\n");
    printf("Input  : mindexx.dat (if model='pcl')\n");
    printf("       : mindexpx.dat, pccfpx.dat (if model='pcx')\n");
    printf("       : mindexp.dat, pccf_all.dat (if model='pc')\n");
    printf("       : mindexp.*.dat, pccfp.*.dat (if model='pcs')\n");
    printf("Output : chain.dat, pchain.dat, mapparam.dat, datavars.dat, parampccfs.dat\n");
    printf("         pmeans.dat, pvars.dat, fmeans.dat, fvars.dat\n");
    printf("       : pdens.dat and pdens_log.dat (if -q)\n");
    printf("================================================================================\n");
    exit(0);
    return 0;
}

///  Main program: Bayesian inference of a few standard function types
int main (int argc, char *argv[])
{
    /// Set the defaults, where necessary
    const char* modeltype=MODELTYPE;
    const char* liktype=LIKTYPE;
    double likParam_dbl=LIKPARAM_DBL;
    int likParam_int=LIKPARAM_INT;
    const char* priortype=PRIORTYPE;
    double priora=PRIORA;
    double priorb=PRIORB;
    const char* datanoise_input;
    const char* chaininitfile;
    const char* xfile=XFILE;
    const char* yfile=YFILE;
    char* xgridfile;
    int pdim=PDIM;
    int order=ORDER;
    const char* rndindfile;
    const char* pdftype=PDFTYPE;
    const char* pctype=PCTYPE;
    double datanoise=DATANOISE;
    int nmcmc=NMCMC;
    double mcmcgamma=MCMCGAMMA;
    char* pgridfile;
    char* pchainfile;
    int fburn=FBURN;
    int nstep=NSTEP;
    char* fixindnomfile;

    bool tflag=false;
    bool eflag=false;
    bool rflag=false;
    bool qflag=false;
    bool pflag=false;
    bool jflag=false;
    bool zflag=false;
    bool vflag=false;

    /// Parse input arguments
    int cc;
    while ((cc=getopt(argc,(char **)argv,"hf:l:w:k:i:a:b:j:x:y:t:d:o:r:s:c:e:m:g:q:p:u:n:zv:"))!=-1){
        switch (cc) {
            case 'h':
                usage();
                break;
            case 'f':
                modeltype = optarg;
                break;
            case 'l':
                liktype = optarg;
                break;
            case 'w':
                likParam_dbl = strtod(optarg, (char **)NULL);
                break;
            case 'k':
                likParam_int = strtol(optarg, (char **)NULL,0);
                break;
            case 'i':
                priortype=optarg;
                break;
            case 'a':
                priora = strtod(optarg, (char **)NULL);
                break;
            case 'b':
                priorb = strtod(optarg, (char **)NULL);
                break;
            case 'j':
                chaininitfile = optarg;
                jflag=true;
                break;
            case 'x':
                xfile = optarg;
                break;
            case 'y':
                yfile = optarg;
                break;
            case 't':
                xgridfile = optarg;
                tflag=true;
                break;
            case 'd':
                pdim =  strtol(optarg, (char **)NULL,0);
                break;
            case 'o':
                order =  strtol(optarg, (char **)NULL,0);
                break;
            case 'r':
                rndindfile = optarg;
                rflag=true;
                break;
            case 's':
                pdftype=optarg;
                break;
            case 'c':
                pctype=optarg;
                break;
            case 'e':
                datanoise_input=optarg;
                eflag=true;
                break;
            case 'm':
                nmcmc =  strtol(optarg, (char **)NULL,0);
                break;
            case 'g':
                mcmcgamma = strtod(optarg, (char **)NULL);
                break;
            case 'q':
                pgridfile=optarg;
                qflag=true;
                break;
            case 'p':
                pchainfile = optarg;
                pflag=true;
                break;
            case 'u':
                fburn =  strtol(optarg, (char **)NULL,0);
                break;
            case 'n':
                nstep =  strtol(optarg, (char **)NULL,0);
                break;
            case 'z':
                zflag = true;
                break;
            case 'v':
                fixindnomfile=optarg;
                vflag=true;
                break;
            default :
                break;
        }
    }



    /// Read datafiles
    Array2D<double> xdata;
    read_datafileVS(xdata,xfile);
    Array1D<Array1D<double> > ydata;
    read_datafileVS(ydata,yfile);
    int nx=xdata.XSize();

    // Define arrays
    Array2D<double> xgrid,pgrid,pchain,chaininit;


    // Data noise indicator
    int dataNoiseInference;
    Array1D<double> datanoise_array;


    /*******************************************************************************************************/
    // Dump the input information
    printf("==============================================================\n") ;
    printf("model_inf() settings:\n");
    printf("---------------------\n");
    printf("Forward model name                     : %s \n",modeltype);
    printf("Likelihood type                        : %s \n",liktype);
    printf("Prior type                             : %s \n",priortype);
    printf("Xdata file                             : %s \n",xfile);
    printf("Ydata file                             : %s \n",yfile);
    printf("Model parameter dim                    : %d \n",pdim);
    printf("Model parameter PC order               : %d \n",order);
    printf("Model parameter PDF type               : %s \n",pdftype);
    printf("Model parameter PC type                : %s \n",pctype);

    if (jflag){
         read_datafileVS(chaininit,chaininitfile);
    }
    else{
        chaininit.Resize(0,2);
    }

    if (rflag){
        printf("Indices of randomized model parameters : %s \n",rndindfile);
    }

    if (tflag){
        read_datafileVS(xgrid,xgridfile);
        printf("Xgrid file for predictions             : %s\n", xgridfile) ;
    }
    else{
        read_datafileVS(xgrid,xfile);
        printf("Xgrid file for predictions             : %s\n", xfile) ;
    }

    if (qflag){
        printf("Posterior will be evaluated at         : %s\n", pgridfile) ;
        read_datafileVS(pgrid,pgridfile);
    }

    if (pflag){
        printf("No MCMC; only postprocessing from file : %s \n", pchainfile) ;
        read_datafileVS(pchain,pchainfile);
    }
    else{
        printf("Number of MCMC samples                 : %d \n",abs(nmcmc));
        if (nmcmc>0)
            printf("Gamma parameter for adaptive MCMC               : %lg \n",mcmcgamma);
        else if (nmcmc<0)
            printf("Non-adaptive MCMC is requested.\n");

    }

    if (eflag){
        if (isalpha(datanoise_input[0])) {
            dataNoiseInference=0;
            read_datafileVS(datanoise_array,datanoise_input);
            printf("Data noise is fixed and read from file %s. \n", datanoise_input) ;

            assert(nx==datanoise_array.Length());
        }
        else
        {

            double datanoise = strtod(datanoise_input, (char **)NULL);

            if (datanoise<0.0){
                dataNoiseInference=1;
                printf("Data noise will be inferred. \n") ;
                datanoise_array.Resize(nx,-datanoise);
            }
            else{
                dataNoiseInference=0;
                printf("Data noise will be fixed at %lg. \n", datanoise) ;
                datanoise_array.Resize(nx,datanoise);
            }
        }

    }
    else{
        dataNoiseInference=2;
        printf("Logarithm of data noise will be inferred. \n") ;
        datanoise_array.Resize(nx,log(datanoise));

    }

    /*******************************************************************************************************/


    // Get xgrid and nburn sizes
    int nxgr=xgrid.XSize();
    int nburn=(int) abs(nmcmc)/fburn;

    /// Read the indices of randomized parameters
    Array1D<int> rndInd;
    if (rflag){
        Array2D<int> rndind2d;
        read_datafileVS(rndind2d,rndindfile);
        assert(rndind2d.YSize()==1);
        array2Dto1D(rndind2d,rndInd);
    }
    else{
//        rndInd.Resize(0);
        for (int i=0;i<pdim;i++)
            rndInd.PushBack(i);
    }

    // Read the indices of fixed parameters
    Array2D<double> fixIndNom(0,2);
    if (vflag){
        read_datafileVS(fixIndNom,fixindnomfile);
        assert (fixIndNom.YSize()==2);
    }

    // Output containers
    Array1D<double> mapparam,pmean_map,pvar_map, fmean_map,fvar_map;
    Array1D<double> datavar_map;
    Array1D<double> p_postave_mean(pdim), p_postave_var(pdim), p_postvar_mean(pdim);
    Array2D<double> f_postsam_mean(nxgr,0);
    Array1D<double> f_postave_mean(nxgr), f_postave_var(nxgr), f_postvar_mean(nxgr);
    Array1D<double> postave_datavar;
    Array2D<double> pmeans,pvars,fmeans,fvars,datavars,paramPCcfs;


    map<string, Array2D<double> (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *) > func_dict;
    func_dict["prop"]            = Func_Prop;
    func_dict["prop_quad"]       = Func_PropQuad;
    func_dict["exp"]             = Func_Exp;
    func_dict["exp_quad"]        = Func_ExpQuad;
    func_dict["const"]           = Func_Const;
    func_dict["linear"]          = Func_Linear;
    func_dict["bb"]              = Func_BB;
    func_dict["heat_transfer1"]  = Func_HT1;
    func_dict["heat_transfer2"]  = Func_HT2;
    func_dict["frac_power"]      = Func_FracPower;
    func_dict["exp_sketch"]      = Func_ExpSketch;
    func_dict["inp"]             = Func_Inputs;
    func_dict["pcl"]             = Func_PCl;
    func_dict["pcx"]             = Func_PCx;
    func_dict["pc"]              = Func_PC;
    func_dict["pcs"]             = Func_PCs;


    if (func_dict.count(modeltype)==0){
        cout << "Model type " << modeltype << " is not found. Exiting." << endl;
        exit(1);
    }
    Array2D<double> (*forwardFunc)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *);
    forwardFunc=func_dict[modeltype];

    int nf=1;
    Array1D< Array2D<double> (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&, void *) > forwardFuncs(nf,NULL);
    forwardFuncs(0)=forwardFunc;
    //forwardFuncs(1)=func_dict["prop"];

    Array1D<double> chstart,chsig;
    getCol(chaininit, 0, chstart);
    getCol(chaininit, 1, chsig);
    bool optimflag=zflag;

    //int seed=SEED;
    srand (time(NULL));
    int seed=rand();
    cout << "Random seed : " << seed << endl;

    int pred_mode=0;
    void* funcinfo=(void*) &pred_mode;
    /// Run the inference
    infer_model(forwardFuncs,funcinfo,
        liktype,
        priortype,priora,priorb,
        xdata,ydata, xgrid,
        dataNoiseInference,datanoise_array,
        pdim,order,rndInd,fixIndNom,pdftype,pctype,
        seed,nmcmc,mcmcgamma, optimflag,chstart, chsig,
        likParam_dbl,likParam_int,
        pgrid, pchain, nburn, nstep,
        mapparam, datavar_map,
        pmean_map,pvar_map,
        fmean_map,fvar_map,
        postave_datavar,
        p_postave_mean,p_postave_var,p_postvar_mean,
        f_postsam_mean,f_postave_mean,f_postave_var,f_postvar_mean,
        paramPCcfs);

    /// Write the outputs
    if (!pflag){
        write_datafile(pchain,"pchain.dat");
        write_datafile_1d(mapparam,"mapparam.dat");
    }

    array1Dto2D(p_postave_mean, pmeans);
    pmeans.insertCol(pmean_map,1);
    write_datafile(pmeans,"pmeans.dat");
    array1Dto2D(p_postave_var, pvars);
    pvars.insertCol(p_postvar_mean,1);
    pvars.insertCol(pvar_map,2);
    write_datafile(pvars,"pvars.dat");

    array1Dto2D(f_postave_mean, fmeans);
    fmeans.insertCol(fmean_map,1);
    write_datafile(fmeans,"fmeans.dat");
    array1Dto2D(f_postave_var, fvars);
    fvars.insertCol(f_postvar_mean,1);
    fvars.insertCol(fvar_map,2);
    write_datafile(fvars,"fvars.dat");
    array1Dto2D(postave_datavar,datavars);
    datavars.insertCol(datavar_map,1);
    write_datafile(datavars,"datavars.dat");

    write_datafile(f_postsam_mean,"fmeans_sams.dat");
    write_datafile(paramPCcfs,"parampccfs.dat");


  return 0;
}

