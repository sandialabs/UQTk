#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.3
#                          Copyright (2023) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
#     retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is open source software: you can redistribute it and/or modify
#     it under the terms of BSD 3-Clause License
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     BSD 3 Clause License for more details.
#
#     You should have received a copy of the BSD 3 Clause License
#     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
#
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
#=====================================================================================

import argparse
import os
import sys
import numpy as np
import math
if (sys.version_info.major==2):
    import cPickle as pick
elif (sys.version_info.major==3):
    import pickle as pick
else:
    print("Only Python 2 or 3 are supported. Exiting.")
    sys.exit()


from model import model

uqtkbin=os.environ['UQTK_INS']+"/bin/"

#############################################################
#############################################################
#############################################################

def model_pc(modelParam, pcparams):
    """PC surrogate evaluator"""

    #print "Running the surrogate model with parameters ", pcparams

    np.savetxt('mindex.dat',pcparams[0],fmt='%d')
    np.savetxt('pccf.dat',pcparams[1])
    pctype=pcparams[2]

    np.savetxt('xdata.dat',modelParam)
    cmd="pce_eval -x'PC_mi' -f'pccf.dat' -s"+pctype+" -r'mindex.dat' > fev.log"
    print("Running %s" % cmd)
    os.system(uqtkbin+cmd)
    pcoutput=np.loadtxt('ydata.dat')
    return pcoutput

#######################################################################################
#######################################################################################
#######################################################################################

def main(argv):

    ## Parse input arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outs', type=int,nargs='*',help="Range of indices of requested outputs (count from 1)")

    parser.add_argument("-r", "--regime",   dest="run_regime",   type=str,   default='online_example', help="Run regime",  choices=['online_example', 'online_bb', 'offline_prep', 'offline_post'])
    parser.add_argument("-p", "--pdom",     dest="pdomain_file", type=str,   default=None, help="Parameter domain file")
    parser.add_argument("-c", "--pcfile",   dest="inpc_file",    type=str,   default=None, help="Input PC coef. file")
    parser.add_argument("-d", "--pcdim",    dest="in_pcdim",     type=int,   default=None, help="Input PC stoch. dimension")
    parser.add_argument("-x", "--pctype",   dest="pc_type",      type=str,   default='LU',    help="PC type",     choices=['HG','LU','LU_N','LG','JB','SW'])
    parser.add_argument("-o", "--pcord",    dest="in_pcord",     type=int,   default=1,       help="Input PC order")
    parser.add_argument("-m", "--method",   dest="fit_method",   type=str,   default='proj',  help="Surrogate construction method", choices=['proj','lsq','bcs'])
    parser.add_argument("-s", "--sampl",    dest="sam_method",   type=str,   default='quad',  help="Sampling method", choices=['quad','rand'])
    parser.add_argument("-n", "--nqd",      dest="num_pt",       type=int,   default=7,       help="Number of quadrature points per dim (if sampl=quad, sparsity=full), \
                                                                                                    Level of quadrature (if sampl=quad, sparsity=sparse), \
                                                                                                    Number of training points (if sampl=rand)")
    parser.add_argument("-v", "--nval",     dest="num_val",      type=int,   default=0,       help="Number of validation (testing) points; can be 0") #  (default: %(default)s)
    parser.add_argument("-f", "--sparsity", dest="sp_type",      type=str,   default='full',  help="Sparsity type", choices=['full', 'sparse'])
    parser.add_argument("-t", "--outord",   dest="out_pcord",    type=int,   default=3,       help="Output PC order")
    parser.add_argument("-i", "--pred",     dest="pred_mode",    type=str,   default='ms',    help="Prediction mode (whether to compute mean only, mean+stdev, or mean+stdev+cov)", choices=['m','ms','msc'])
    parser.add_argument("-e", "--tol",      dest="tolerance",    type=float, default=1.e-3,   help="Tolerance (currently for method=bcs only)")
    parser.add_argument("-z", "--cleanup",  dest="cleanup", default=False, action="store_true" , help="Flag to cleanup after (be careful: removes *log and *dat files)")
    parser.add_argument("-y", "--normalize",  dest="normalize", default=False, action="store_true" , help="Normalize the model outputs before constructing surrogates")
    args = parser.parse_args()

    ## Flags for input checks
    pflag=False
    cflag=False
    dflag=False


    ## Hardwired names
    input_train='ptrain.dat'
    input_val='pval.dat'
    qd_train='qtrain.dat'
    wg_train='wtrain.dat'
    qd_val='qval.dat'
    output_train='ytrain.dat'
    output_val='yval.dat'


    ## Argument compatibility checks
    if args.pdomain_file!=None:
        pflag=True
    if args.inpc_file!=None:
        cflag=True
    if args.in_pcdim!=None:
        dflag=True

    ## Sanity check
    if (args.run_regime!="offline_prep" and args.fit_method=='proj' and args.sam_method!='quad'):
        print("Projection method requires quadrature sampling. Exiting.")
        sys.exit()

    ## Arguments that may be overwritten
    in_pcord = args.in_pcord
    pc_type  = args.pc_type

    ## Organize input PC structure
    if (pflag and cflag):
        print("Need to provide input domain or input PC coef file, not both. Exiting.")
        sys.exit()
    elif (pflag and dflag):
        print("Need to provide input domain or input dimensionality, not both. Exiting.")
        sys.exit()
    elif (int(pflag)+int(cflag)+int(dflag)==0):
        in_pcdim=2 # default if no argument is given
        pcf_all=np.vstack((np.zeros((in_pcdim,)),np.eye(in_pcdim)))
    elif pflag:
        print("Input order (-o) and pctype (-x) will be overwritten, if domain file %s is given. " % args.pdomain_file)
        # Load parameter domain file
        if (os.path.isfile(args.pdomain_file)):
            pdom=np.loadtxt(args.pdomain_file).reshape(-1,2)
            assert(pdom.shape[1]==2)
            pcf_all=np.vstack((0.5*(pdom[:,1]+pdom[:,0]),np.diag(0.5*(pdom[:,1]-pdom[:,0]))))
            in_pcdim=pdom.shape[0]
            in_pcord=1
            if pc_type[:2] != "LU":
                pc_type="LU"

            for i in range(pdom.shape[0]):
                if(pdom[i,0]>pdom[i,1]):
                    print("Error: The domain file %s contains wrong bounds. Check the row number %d. Exiting." % (args.pdomain_file,i+1))
                    sys.exit()
        else:
            print("Error: The requested domain file %s does not exist. Exiting." % args.pdomain_file)
            sys.exit()
    elif (cflag):
        if (os.path.isfile(args.inpc_file)):
            pcf_all=np.atleast_2d(np.loadtxt(args.inpc_file))
            in_pcdim=args.in_pcdim
            # TODO sanity checks on in_pcdim and in_pcord with the size of pcf_all
        else:
            print("Error: The requested input PC coefficient file %s does not exist. Exiting." % args.inpc_file)
            sys.exit()
    elif (dflag):
        in_pcdim=args.in_pcdim
        pcf_all=np.vstack((np.zeros((in_pcdim,)),np.eye(in_pcdim)))
    else:
        print("If this message appears, there must be a bug in the code. Exiting.")
        sys.exit()

    ## Sanity check to ensure projection uses quadrature points
    if args.run_regime!="offline_prep" and args.fit_method=="proj" and args.sam_method!="quad":
        print("Projection requires quadrature sampling. Exiting.")
        sys.exit(1)

    ## Get the dimensions
    npar = pcf_all.shape[1]
    npc  = pcf_all.shape[0] # not used

    ## Print the inputs for reference
    print("Run regime                        %s" % args.run_regime)
    if (cflag):
        print("Input PC coefficient file         %s" % args.inpc_file)
    if (pflag):
        print("Input parameter domain file       %s" % args.pdomain_file)
    print("Input PC dim                      %d" % in_pcdim)
    print("Input PC order                    %d" % in_pcord)
    print("PC type                           %s" % pc_type)
    print("The number of input parameters    %d" % npar)
    if args.run_regime!="offline_prep":
        print("Method                            %s" % args.fit_method)
    print("Sampling method                   %s" % args.sam_method)
    print(" with parameter                   %d" % args.num_pt)
    print("Number of validation points       %d" % args.num_val)
    if (args.sam_method=="quad"):
        print("Sparsity type                     %s" % args.sp_type)
    if args.run_regime!="offline_prep":
        print("Output PC order                   %d" % args.out_pcord)

    ## (1) Generate sample points for online or offline_prep regimes
    if args.run_regime!="offline_post":
        print("#####################################################################")
        print("######################## Generate input samples #####################")
        print("#####################################################################")
        print("#### Generating training samples")

        if args.sam_method=="quad":
            cmd="generate_quad -d"+str(in_pcdim)+"  -g"+pc_type+" -x"+args.sp_type+" -p"+str(args.num_pt)+" > gq.log"
            print("Running %s" % cmd)
            os.system(uqtkbin+cmd)
            inqdp=np.loadtxt('qdpts.dat')
            inqdw=np.loadtxt('wghts.dat')
            np.savetxt(wg_train,inqdw)


        elif args.sam_method=="rand":
            seedval=np.random.randint(100, size=10)

            cmd="pce_rv -w PCvar -d"+str(in_pcdim)+"  -p"+str(in_pcdim)+" -x"+pc_type+" -n"+str(args.num_pt)+" -s "+str(seedval[0])+" > pcrv.log"
            print("Running %s" % cmd)
            os.system(uqtkbin+cmd)

            inqdp=np.loadtxt('rvar.dat')

        else:
            print("Error: Sampling method is not recognized. Should be 'quad' or 'rand'. Exiting.")
            sys.exit()

        np.savetxt(qd_train,inqdp)
        npt=inqdp.shape[0]
        print("Germ samples for training are in %s in a format %d x %d " % (qd_train,npt,in_pcdim))

        ## Evaluate input PCs at quadrature points
        np.savetxt('xdata.dat',inqdp)
        inpar=np.empty((npt,npar))
        for i in range(npar):
            np.savetxt('pcf.dat',pcf_all[:,i])
            cmd="pce_eval -x PC -s"+pc_type+" -o"+str(in_pcord)+" -f 'pcf.dat' > pcev.log"
            print("Running %s" % cmd)
            os.system(uqtkbin+cmd)
            inpar[:,i]=np.loadtxt('ydata.dat')

        np.savetxt(input_train,inpar)
        print("Parameter samples for training are in %s in a format %d x %d " % (input_train,npt,npar))

        ## Generate points, if requested, for the validation of the surrogate
        seedval=np.random.randint(100, size=10)
        if args.num_val>0:
            print("#### Generating validation samples")
            cmd="pce_rv -w PCvar -d "+str(in_pcdim)+ " -n "+str(args.num_val)+ " -p "+str(in_pcdim)+" -x "+pc_type+" -s "+str(seedval[1])+" > pcrv.log"
            print("Running %s" % cmd)
            os.system(uqtkbin+cmd)
            qpar_val=np.loadtxt('rvar.dat')
            np.savetxt(qd_val,qpar_val)
            print("Germ samples for validation are in %s in a format %d x %d " % (qd_val,args.num_val,in_pcdim))


            inpar_val=np.empty((args.num_val,npar))
            for i in range(npar):
                np.savetxt('pcf.dat',pcf_all[:,i])
                cmd="cp rvar.dat xdata.dat"
                os.system(cmd)
                cmd="pce_eval -x PC -s"+pc_type+" -o"+str(in_pcord)+" -f 'pcf.dat' > pcev.log"
                print("Running %s" % cmd)
                os.system(uqtkbin+cmd)
                inpar_val[:,i]=np.loadtxt('ydata.dat')


            np.savetxt(input_val,inpar_val)
            print("Parameter samples for validation are in %s in a format %d x %d " % (input_val,args.num_val,npar))

        ## Exit if only sample preparation is required
        if args.run_regime=="offline_prep":
            print("Preparation of samples is done.")
            sys.exit()

    print("#####################################################################")
    print("######################## Load input samples #########################")
    print("#####################################################################")
    ## (2) Load sample points for online or offline_post regimes
    ptrain=np.loadtxt(input_train).reshape(-1,npar)
    inqdp=np.loadtxt(qd_train).reshape(-1,in_pcdim)
    if (args.num_val>0):
        pval=np.loadtxt(input_val).reshape(-1,npar)
        qpar_val=np.loadtxt(qd_val).reshape(-1,in_pcdim)

    npt=ptrain.shape[0]
    if (args.num_val>0):
        args.num_val=pval.shape[0]

    print("Number of training points for surrogate construction   : %d" % npt)
    print("Number of validation points for surrogate construction : %d" % args.num_val)

    print("#####################################################################")
    print("######################## Evaluate/Load forward model  ###############")
    print("#####################################################################")
    ## (3) Get model outputs

    # Run the model online or....
    if args.run_regime=="online_example":
        ytrain=model(ptrain)
        if (args.num_val>0):
            yval=model(pval)

    elif args.run_regime=="online_bb":
        # TODO check model.x existence
        os.system('./model.x '+input_train+' '+output_train)
        ytrain=np.loadtxt(output_train).reshape(npt,-1)
        if (args.num_val>0):
            os.system('./model.x '+input_val+' '+output_val)
            yval=np.loadtxt(output_val).reshape(args.num_val,-1)

    # ...or read the results from offline simulations
    elif args.run_regime=="offline_post":
        ytrain=np.loadtxt(output_train).reshape(npt,-1)
        if (args.num_val>0):
            yval=np.loadtxt(output_val).reshape(args.num_val,-1)



    # Read the number of output observables or the number of values of deisgn parameters (e.g. location, time etc..)
    nout_all=ytrain.shape[1]
    print("Number of output observables of the model : %d" % nout_all)

    # Normalize ytrain and yval if asked
    if args.normalize:
        yscale = np.max(np.abs(ytrain), axis=0)
    else:
        yscale = np.ones(nout_all)

    ytrain /= yscale
    if (args.num_val > 0):
        yval /= yscale

    print("#####################################################################")
    print("######################## Construct PC surrogates  ###################")
    print("#####################################################################")
    ## (4) Obtain the PC surrogate using model simulations
    if len(args.outs)==0:
        outrange=range(nout_all)
    elif len(args.outs)==1:
        outrange=range(args.outs[0]-1,args.outs[0])
    elif len(args.outs)==2:
        outrange=range(min(args.outs)-1,max(args.outs))
    else:
        print("The number of free arguments can not be more than two. Exiting.")
        sys.exit()

    ## TODO put sanity checks in
    nout = len(outrange)
    print("Number of output observables being analyzed : %d" % nout)


    ## Empty arrays and lists to store results
    pccf_all=[]
    mindex_all=[]
    varfrac_all=[]
    ccov_all=[]
    allsens_main=np.empty((nout,in_pcdim))
    allsens_total=np.empty((nout,in_pcdim))
    allsens_joint=np.empty((nout,in_pcdim, in_pcdim))
    ytrain_pc=np.empty((npt,nout))
    yval_pc=np.empty((args.num_val,nout))
    errcheck_pc=np.empty((npt,nout))
    errcheck_val_pc=np.empty((args.num_val,nout))

    err_training=np.empty((nout,))
    err_val=np.empty((nout,))

    ## Generate PC multiindex
    cmd="gen_mi -x'TO' -p"+str(args.out_pcord)+" -q"+str(in_pcdim)+" > gmi.log; mv mindex.dat mi.dat"
    print("Running %s" % cmd)
    os.system(uqtkbin+cmd)
    #os.system('cp custom_mindex.dat mi.dat')
    mi=np.loadtxt('mi.dat')

    npc=mi.shape[0]

    inpar=np.loadtxt(input_train)
    xcheck=inqdp.copy()
    if (args.num_val>0):
        inpar_val=np.loadtxt(input_val)
        qpar_val=np.loadtxt(qd_val).reshape(args.num_val,-1)
        xcheck=np.vstack((inqdp,qpar_val))
    ncheck=xcheck.shape[0]

    # Loop over all output observables/locations
    i = 0
    for j in outrange:
        ################################
        # (4a) Build PC surrogate
        print("##################################################")
        print("Building PC for observable %d / %d" % (j+1,nout_all))
        np.savetxt('ydata.dat',ytrain[:,j])
        inqdp=np.loadtxt(qd_train)

        if args.fit_method=="proj":

            inqdw=np.loadtxt(wg_train)
            np.savetxt('qdpts.dat',inqdp)
            np.savetxt('wghts.dat',inqdw)

            cmd="pce_resp -x"+pc_type+" -o"+str(args.out_pcord)+" -d"+str(in_pcdim)+" -e > pcr.log"
            print("Running %s" % cmd)
            os.system(uqtkbin+cmd)

            # Get the PC coefficients and multiindex
            pccf=np.loadtxt('PCcoeff_quad.dat').reshape(-1,)
            mindex=np.loadtxt('mindex.dat',dtype='int').reshape(-1,in_pcdim)

            ycheck_var=np.zeros((ncheck,))
            sig2=0.0
            ccov=np.zeros((mindex.shape[0],mindex.shape[0]))
            erb=np.zeros((ncheck,))

        elif args.fit_method=="bcs":
            bcs_tol=args.tolerance
            np.savetxt('xdata.dat',inqdp)
            np.savetxt('xcheck.dat',xcheck)
            np.savetxt('regparams.dat',np.ones((npc,)))
            cmd='regression -x xdata.dat -y ydata.dat -b PC_MI -s '+pc_type+' -p mi.dat -w regparams.dat -m '+args.pred_mode+' -r wbcs -t xcheck.dat -c '+str(bcs_tol)+' > regr.log'
            print("Running %s" % cmd)
            os.system(uqtkbin+cmd)

            # Get the PC coefficients and multiindex and the predictive errorbars
            pccf=np.loadtxt('coeff.dat').reshape(-1,)
            mindex=np.loadtxt('mindex_new.dat',dtype='int').reshape(-1,in_pcdim)

            ycheck_var=np.zeros((ncheck,))
            sig2=0.0
            ccov=np.zeros((mindex.shape[0],mindex.shape[0]))
            if (args.pred_mode!='m'):
                sig2=max(0.0,np.loadtxt('sigma2.dat'))
                ccov=np.loadtxt('Sig.dat')
                #if os.path.isfile('ycheck_var.dat'):
                ycheck_var=np.loadtxt('ycheck_var.dat')

            erb=np.sqrt(ycheck_var+sig2)


        elif args.fit_method=="lsq":
            np.savetxt('xdata.dat',inqdp)
            np.savetxt('xcheck.dat',xcheck)
            np.savetxt('regparams.dat',np.ones((npc,)))
            cmd='regression -l0 -x xdata.dat -y ydata.dat -b PC_MI -s '+pc_type+' -p mi.dat -m '+args.pred_mode+' -r lsq -t xcheck.dat > regr.log'
            print("Running %s" % cmd)
            os.system(uqtkbin+cmd)

            # Get the PC coefficients and multiindex and the predictive errorbars
            pccf=np.loadtxt('coeff.dat').reshape(-1,)
            mindex=np.loadtxt('mi.dat',dtype='int').reshape(-1,in_pcdim)

            ycheck_var=np.zeros((ncheck,))
            sig2=0.0
            ccov=np.zeros((mindex.shape[0],mindex.shape[0]))
            if (args.pred_mode!='m'):
                sig2=max(0.0,np.loadtxt('sigma2.dat'))
                ccov=np.loadtxt('Sig.dat')
                # This is a cheap hack to tolerate Cholesky error during LSQ when simpler polynomials are given much more complex mindices
                #if os.path.isfile('ycheck_var.dat'):
                ycheck_var=np.loadtxt('ycheck_var.dat')

            erb=np.sqrt(ycheck_var+sig2)


        else:
            print("Method not recognized. Exiting.")
            sys.exit()


        # Scale
        pccf *= yscale[j]
        ccov *= yscale[j]**2
        ytrain[:, j] *= yscale[j]
        erb *= yscale[j]
        if (args.num_val>0):
            yval[:, j] *= yscale[j]

        # Append the results
        pccf_all.append(pccf)
        mindex_all.append(mindex)
        ccov_all.append(ccov)

        ################################

        # (4b) Evaluate the PC surrogate at training and validation points
        print("Evaluating surrogate at %d training points" % npt)
        ytrain_pc[:,i]=model_pc(inqdp,[mindex,pccf,pc_type])
        normf=np.linalg.norm(ytrain[:,j])
        if normf!=0.0:
            err_training[i]=np.linalg.norm(ytrain[:,j]-ytrain_pc[:,i])/normf
        else:
            err_training[i]=0.0
        print("Surrogate relative error at training points : %.12g" % err_training[i])

        errcheck_pc[:,i]=erb[:npt]




        if (args.num_val>0):

            print("Evaluating surrogate at %d validation points" % args.num_val)
            yval_pc[:,i]=model_pc(qpar_val,[mindex,pccf,pc_type])
            normf=np.linalg.norm(yval[:,j])
            if normf!=0.0:
                err_val[i]=np.linalg.norm(yval[:,j]-yval_pc[:,i])/normf
            else:
                err_val[i]=0.0
            print("Surrogate relative error at validation points : %.12g" % err_val[i])
            #np.savetxt('yval_pc.'+str(i+1)+'.dat',yval_pc)
            errcheck_val_pc[:,i]=erb[npt:]



        ################################

        # (4c) Compute sensitivities
        np.savetxt('PCcoeff.dat',pccf)
        cmd="pce_sens -m'mindex.dat' -f'PCcoeff.dat' -x"+pc_type+" > pcsens.log"
        print("Running %s" % cmd)
        os.system(uqtkbin+cmd)
        allsens_main[i,:]=np.loadtxt('mainsens.dat')
        allsens_total[i,:]=np.loadtxt('totsens.dat')
        allsens_joint[i,:,:]=np.loadtxt('jointsens.dat')
        print("Sum of main sensitivities  : %.12g" % sum(allsens_main[i,:]))
        print("Sum of total sensitivities : %.12g" % sum(allsens_total[i,:]))
        varfrac=np.loadtxt('varfrac.dat')
        varfrac_all.append(varfrac)

        i+=1
    ############################################################################


    ## Results container
    if(args.num_val>0):
        results = {'outs':(outrange),'training':(inqdp,inpar,ytrain,ytrain_pc,errcheck_pc),'validation':(qpar_val,inpar_val,yval,yval_pc,errcheck_val_pc),'pcmi':(pccf_all,mindex_all,varfrac_all,ccov_all,args.pc_type),'sens':(allsens_main,allsens_total,allsens_joint),'err':(err_training,err_val)}
    else:
        results = {'outs':(outrange),'training':(inqdp,inpar,ytrain,ytrain_pc,errcheck_pc),'pcmi':(pccf_all,mindex_all,varfrac_all,ccov_all,args.pc_type),'sens':(allsens_main,allsens_total,allsens_joint),'err':(err_training,)}

    ## Save results
    pick.dump(results,open('results.pk','wb'),-1)


    # Cleanup of unneeded leftovers
    if args.cleanup:
        del_cmd='rm -rf *dat *log'
        os.system(del_cmd)

if __name__ == "__main__":
   main(sys.argv[1:])
