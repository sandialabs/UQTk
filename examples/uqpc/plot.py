#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.2
#                          Copyright (2022) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

from __future__ import print_function
import argparse
import os
import sys
try:
    import numpy as np
except ImportError:
    print("Numpy was not found.")

try:
    import matplotlib
except ImportError:
    print("Matplotlib was not found.")
import math



sys.path.append(os.environ['UQTK_INS'])
import PyUQTk.plotting.surrogate as ss
import PyUQTk.plotting.inout as io
import PyUQTk.utils.func as uf

from pylab import *

if (sys.version_info.major==2):
    import cPickle as pick
elif (sys.version_info.major==3):
    import pickle as pick

rc('legend',loc='best', fontsize=22)
rc('lines', linewidth=2, color='r')
rc('axes',linewidth=3,grid=True,labelsize=28)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)


#############################################################
#############################################################

## Parsing the inputs
usage_str='           Reads results.pk and makes visualizations.\n\
           First argument defines the plot type. \n\
           Try "plot.py <plot_type> -h" for help on additional arguments under each plot_type. \n\
           Parameter names and output names can be read in pnames.txt and outnames.txt (if files not present, will use generic names).\n\
           Many options are quite experimental yet, and not optimal visually. \n\
           One can use this script as an example of how to unroll results.pk for subsequent plotting on one\'s own.'

parser = argparse.ArgumentParser(description=usage_str,formatter_class=argparse.RawTextHelpFormatter)
arg1_choices=['sens','senscirc','sensmat','dm','idm','1d','2d','mindex','micf','pdf','senserb1','senserb2']
#parser.add_argument('plot_type', type=str,nargs=1,help="Plot type", choices=arg1_choices)
subparsers=parser.add_subparsers()
for arg1 in arg1_choices:
    if arg1=='sens':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots sensitivity barplot for all outputs', \
            epilog='Examples: \n"plot.py sens main", \n"plot.py sens total"')
        sbparser.add_argument('sensmode', type=str,nargs='?',help="Sensitivity type",choices=['main','total'])
    elif arg1=='senscirc':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots sensitivity circular plots for all outputs, and averaged as well.', \
            epilog='Example: \n"plot.py senscirc"')
    elif arg1=='sensmat':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots sensitivity matrix for all outputs and for the most important inputs', \
            epilog='Examples: \n"plot.py sensmat main", \n"plot.py sensmat total"')
        sbparser.add_argument('sensmode', type=str,nargs='?',help="Sensitivity type",choices=['main','total'])
    elif arg1=='dm' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots model-vs-surrogate for all outputs', \
            epilog='Examples: \n"plot.py dm training", \n"plot.py dm validation", \n"plot.py dm training validation"')
        sbparser.add_argument('trvals', type=str,nargs='*',help="Data types", choices=['training','validation'])
    elif arg1=='idm' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots model and surrogate versus point id', \
            epilog='Examples: \n"plot.py idm training", \n"plot.py idm validation"')
        sbparser.add_argument('trval', type=str,nargs='?',help="Data type", choices=['training','validation'])
    elif arg1=='1d' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots 1d surrogate (the rest of parameters, if any, at nominal) versus data, for all outputs', \
            epilog='Examples: \n"plot.py 1d 0 training"')
        sbparser.add_argument('d1', type=int,nargs='?',help="First dimension (count from 0)")
        sbparser.add_argument('trval', type=str,nargs='?',help="Data type", choices=['training','validation'])
    elif arg1=='2d' :
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Plots 2d surrogate (the rest of parameters, if any, at nominal) versus data, for all outputs', \
            epilog='Examples: \n"plot.py 2d 0 1 training"')
        sbparser.add_argument('d1', type=int,nargs='?',help="First dimension (count from 0)")
        sbparser.add_argument('d2', type=int,nargs='?',help="Second dimension (count from 0)")
        sbparser.add_argument('trval', type=str,nargs='?',help="Data type", choices=['training','validation'])
    elif arg1=='mindex':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
            description='Visualizes the multiindex for all outputs.', \
            epilog='Example: \n"plot.py mindex"')
    elif arg1=='micf':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
        description='Plots the multiindex for all outputs in a different way, meaningful only for 2d and 3d.', \
        epilog='Example: \n"plot.py micf"')
    elif arg1=='pdf':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
        description='Plots the PDF of the output. Sampling parameter is hardwired.', \
        epilog='Example: \n"plot.py pdf"')
    elif arg1=='senserb1':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
        description='Computes sensitivities with errorbars. Not tested enough. Some hardwired parameters. \n\
                     Requires uq_pc.py method (-m) lsq or bcs and prediction mode (-i) msc.', \
        epilog='Example: \n"plot.py senserb1"')
    elif arg1=='senserb2':
        sbparser = subparsers.add_parser(arg1,formatter_class=argparse.RawTextHelpFormatter,\
        description='Plots the sensitivities with errorbars. Needs to be run after "plot.py senserb1".', \
        epilog='Example: \n"plot.py senserb2"')
    else:
        print("The code should not get here. Please double check the arguments. Exiting.")
        sys.exit()

args = parser.parse_args()
plot_type=sys.argv[1]

## Read surrogate results
results=pick.load(open('results.pk', 'rb'))
print("results.pk dictionary contents : ", results.keys())

## Get basic dimensions
nout=results['training'][2].shape[1]
ndim=results['sens'][0].shape[1]
pctype=results['pcmi'][4]
outrange=results['outs']
print("Dimensionality : %d" % ndim)
print("Num of outputs : %d" % nout)

## Parameter names file, if any. Must have ndim rows
if os.path.exists('pnames.txt'):
    with open('pnames.txt') as f:
        pnames = f.read().splitlines()
        assert(len(pnames)==ndim)
else:
    pnames=['Param # '+str(i) for i in range(1,ndim+1)]

## Output names file, if any. Must have nout rows
if os.path.exists('outnames.txt'):
    with open('outnames.txt') as f:
        outnames = f.read().splitlines()
        assert(len(outnames)==nout)
else:
    outnames=['Output # '+str(i) for i in range(1,nout+1)]


if(plot_type=='sens'):
    allsens_main,allsens_total,allsens_joint=results['sens']

    sensmode=args.sensmode
    if sensmode=='main':
        sensdata=allsens_main
    elif sensmode=='total':
        sensdata=allsens_total

    pars=range(ndim)
    cases=range(nout)

    np.savetxt('allsens_'+sensmode+'.dat',sensdata)


    print('===================================================')
    print("Plotting %s sensitivities across all output QoIs (bar plot)" % sensmode)
    ss.plot_sens(sensdata[outrange],pars,cases,vis="bar",reverse=False,par_labels=pnames,case_labels=outnames,ncol=5,grid_show=False, xlbl='', topsens=10, showplot=False)
    os.system('mv sens.eps sens_'+sensmode+'.eps')
    print('Sensitivities are reported in allsens_%s.dat' % sensmode)
    print('Created file sens_%s.eps' % sensmode)

elif(plot_type=='senscirc'):
    allsens_main,allsens_total,allsens_joint=results['sens']


    npar=allsens_main[0].shape[0]
    allsens_main_ave=np.zeros((npar,))
    allsens_joint_ave=np.zeros((npar,npar))
    print('Plotting main sand joint sensitivities (circular plots)')
    for i in range(len(outrange)):
        iout=outrange[i]
        print('======== Output # ',str(iout+1),' (QoI ',outnames[iout],') =========')
        ss.plot_senscirc(outnames[iout], allsens_main[i],allsens_joint[i],pnames,showplot=False)
        os.system('mv senscirc.eps senscirc_output_'+str(iout+1)+'.eps')
        print('Created file senscirc_output_',str(iout+1),'.eps')

        allsens_main_ave+=allsens_main[i]/npar
        allsens_joint_ave+=allsens_joint[i]/npar

    print('Plotting averaged main sand joint sensitivities (circular plots)')
    ss.plot_senscirc('', allsens_main_ave,allsens_joint_ave,pnames,showplot=False)
    os.system('mv senscirc.eps senscirc_ave.eps')
    print('Created file senscirc_ave.eps')


elif(plot_type=='sensmat'):
    allsens_main,allsens_total,allsens_joint=results['sens']

    sensmode=sys.argv[2]
    if sensmode=='main':
        sensdata=allsens_main
    elif sensmode=='total':
        sensdata=allsens_total


    pars=range(ndim)
    cases=range(nout)

    print('===================================================')
    print("Plotting ",sensmode," sensitivities across all output QoIs (matrix plot)")
    ss.plot_sensmat(sensdata[outrange],pars,cases,vis="bar",reverse=False,par_labels=pnames,case_labels=outnames,showplot=False)
    os.system('mv sensmat.eps sensmat_'+sensmode+'.eps')
    print('Created file sensmat_'+sensmode+'.eps')


elif(plot_type=='dm'):

    # Parse the arguments
    trvals=args.trvals #sys.argv[2:]

    for i in range(len(outrange)):
        iout=outrange[i]
        axes_labels=['Model ('+outnames[iout]+')','Polynomial Surrogate']
        print('======== Output # ',str(iout+1),' (QoI ',outnames[iout],') =========')
        datas=[]
        models=[]
        errorbars=[]
        labels=[]
        for trval in trvals:
            print("Plotting model-vs-surrogate at ",trval," points")
            if trval not in results.keys():
                print(trval, " points are not present in results. Exiting.")
                sys.exit()
            datas.append(results[trval][2][:,iout])
            models.append(results[trval][3][:,i])
            errorbars.append([results[trval][4][:,i],results[trval][4][:,i]])
            labels.append(trval+' points')

        ss.plot_dm(datas,models,errorbars,labels,axes_labels,showplot=False)
        os.system('mv dm.eps dm_output_'+str(iout+1)+'.eps')
        print('Created file dm_output_'+str(iout+1)+'.eps')

elif(plot_type=='idm'):

    # Parse the arguments
    trval=args.trval #sys.argv[2]
    if trval not in results.keys():
        print(trval, " points are not present in results. Exiting.")
        sys.exit()

    for i in range(len(outrange)):
        iout=outrange[i]
        #axes_labels=['Model ('+outnames[iout]+')','Polynomial Surrogate']
        print('======== Output # ',str(iout+1),' (QoI ',outnames[iout],') =========')

        print("Plotting runId-vs-model and runId-vs-surrogate at ",trval," points")
        data=results[trval][2][:,iout]
        model=results[trval][3][:,i]
        errbar=[results[trval][4][:,i],results[trval][4][:,i]]
        #label='trval+' points

        ss.plot_idm(data,model,errbar,sort='none',figname='idm.eps')

        os.system('mv idm.eps idm_output_'+trval+'_'+str(iout+1)+'.eps')
        print('Created file idm_output_'+trval+'_'+str(iout+1)+'.eps')

elif(plot_type=='1d'):
    # Parse the arguments
    d1=args.d1 #int(sys.argv[2])
    trval=args.trval #sys.argv[3]

    if (d1<0 or d1>ndim-1):
        print("The dimension is incorrect. Exiting.")
        sys.exit()
    if trval not in results.keys():
        print(trval, " points are not present in results. Exiting.")
        sys.exit()
    for i in range(len(outrange)):

        iout=outrange[i]
        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        #axes_labels=['Model ('+outnames[iout]+')','Polynomial Surrogate']
        print('======== Output # ',str(iout+1),' (QoI ',outnames[iout],') =========')
        print("Plotting model-vs-surrogate at ",trval," points")
        output=results[trval][2][:,iout]
        input=results[trval][0].reshape(-1,ndim)
        fig=io.plot_xy(input[:,d1],output,pnames[d1], outnames[iout],label=trval+' points',savefig='xy_'+pnames[d1]+'_'+outnames[iout]+'.eps')

        ngr=100
        #prange=np.loadtxt('prange.dat',ndmin=2)
        #input_grid_phys=np.linspace(prange[d,0],prange[d,1],ngr)
        input_grid=np.linspace(-1,1,ngr).reshape(-1,1)

        pc_grid=uf.func(input_grid,'PCmi',[mindex,pcf,pctype])

        fig.gca().plot(input_grid,pc_grid,'r-',ms=11,label='Polynomial fit')
        #xdel=prange[d,1]-prange[d,0]
        #gca().set_xlim(prange[d,0]-0.1*xdel,prange[d,1]+0.1*xdel)
        legend()

        fig.savefig('fit_'+pnames[d1]+'_'+outnames[iout]+'_'+trval+'.eps')
        print('Created file fit_'+pnames[d1]+'_'+outnames[iout]+'_'+trval+'.eps')

elif(plot_type=='2d'):

    # Parse the arguments
    d1=args.d1 #int(sys.argv[2])
    d2=args.d2 #int(sys.argv[3])
    trval=args.trval #sys.argv[4]

    if (d1<0 or d1>ndim-1):
        print("First dimension is incorrect. Exiting.")
        sys.exit()
    if (d2<0 or d2>ndim-1):
        print("Second dimension is incorrect. Exiting.")
        sys.exit()

    if trval not in results.keys():
        print(trval, " points are not present in results. Exiting.")
        sys.exit()

    for i in range(len(outrange)):

        iout=outrange[i]
        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        #axes_labels=['Model ('+outnames[iout]+')','Polynomial Surrogate']
        print('======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') =========')
        print("Plotting model-vs-surrogate at "+trval+" points")
        output=results[trval][2][:,iout]
        input=results[trval][0]
        fig=io.plot_xxy(input[:,d1],input[:,d2],output,pnames, outnames[iout],label=trval+' points',savefig='xxy_'+pnames[d1]+'_'+pnames[d2]+'_'+outnames[iout]+'.eps')


        ngr=100
        #prange=np.loadtxt('prange.dat',ndmin=2)
        #input1_grid_phys=np.linspace(prange[d1,0],prange[d1,1],ngr)
        #input2_grid_phys=np.linspace(prange[d2,0],prange[d2,1],ngr)
        input_grid=np.linspace(-1,1,ngr).reshape(-1,1)
        x,y=np.meshgrid(input_grid,input_grid)
        #xp,yp=np.meshgrid(input1_grid_phys,input2_grid_phys)

        xy=np.vstack([x.ravel(), y.ravel()]).T
        pc_grid=uf.func(xy,'PCmi',[mindex,pcf,pctype]).reshape(ngr,ngr)
        gca().plot_wireframe(x, y, pc_grid, rstride=5,cstride=5,linewidth=0.3)
        #gca().view_init(38,-169)
        fig.savefig('fit_'+pnames[d1]+'_'+pnames[d2]+'_'+outnames[iout]+'_'+trval+'.eps')
        print('Created file fit_'+pnames[d1]+'_'+pnames[d2]+'_'+outnames[iout]+'_'+trval+'.eps')

elif(plot_type=='mindex'):

    print('Plotting multiindices')

    for i in range(len(outrange)):
        iout=outrange[i]
        print('======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') =========')

        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        varfrac=results['pcmi'][2][iout]

        np.savetxt("pcf_output_"+str(iout+1)+".dat",pcf)
        np.savetxt("mindex_output_"+str(iout+1)+".dat",mindex,fmt='%d')
        ss.plot_mindex(mindex[1:,:],varfrac[1:],str(outnames[iout]),showplot=False)
        os.system('mv mindex.eps mindex_output_'+str(iout+1)+'.eps')
        print('Created file mindex_output_'+str(iout+1)+'.eps')

elif(plot_type=='micf'):
    if (ndim>3):
        print("plot_micf utility only visualizes 2d or 3d multiindices. Exiting.")
        sys.exit()

    for i in range(len(outrange)):
        iout=outrange[i]
        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        np.savetxt("pcf_output_"+str(iout+1)+".dat",pcf)
        np.savetxt("mindex_output_"+str(iout+1)+".dat",mindex,fmt='%d')
        ss.plot_micf(mindex,cfs=pcf,showplot=False)
        os.system('mv micf.eps micf_output_'+str(iout+1)+'.eps')
        print('Created file micf_output_'+str(iout+1)+'.eps')

elif(plot_type=='pdf'):
    nsam=100000

    print('Computing and plotting output PDFs')

    for i in range(len(outrange)):
        iout=outrange[i]
        print('======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') =========')

        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        pctype=results['pcmi'][-1]
        order= np.amax(mindex)

        custom_xlabel='Output QoI ('+outnames[iout]+')'

        ss.plot_pcpdf(pctype,mindex,pcf,nsam,custom_xlabel,showplot=False)
        os.system('mv pcdens.eps pcdens_output_'+str(iout+1)+'.eps')
        print('Created file pcdens_output_'+str(iout+1)+'.eps')

elif(plot_type=='senserb1'):
    sensord=1
    senssam=1000
    print('Computing sensitivity errorbars')

    means=np.empty((len(outrange),ndim))
    stds=np.empty((len(outrange),ndim))
    for i in range(len(outrange)):
        iout=outrange[i]
        print('======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') =========')

        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        cfcov=results['pcmi'][3][iout]
        #print results['pcmi'][2]
        #print cfcov, cfcov.shape
        lmat=np.linalg.cholesky(cfcov)
        inputpc=np.vstack((pcf,lmat.T))

        np.savetxt('mi.dat',mindex, fmt='%d')
        np.savetxt('inputpc.dat',inputpc)
        npc=inputpc.shape[1]

        cmd="mv results.pk results_surr.pk"
        os.system(cmd)
        cmd="ln -sf "+os.environ['UQTK_INS']+"/examples/uqpc/model_sens.x model.x"
        os.system(cmd)
        cmd=os.environ['UQTK_INS']+"/examples/uqpc/uq_pc.py -r online_bb -c inputpc.dat -x HG -d "+str(npc)+" -o 1 -m lsq -s rand -f sparse -n "+str(senssam)+" -v 10 -t "+str(sensord)
        os.system(cmd)
        cmd="mv results.pk results_sens_"+str(iout+1)+".pk; mv results_surr.pk results.pk"
        os.system(cmd)


elif(plot_type=='senserb2'):

    print('Plotting sensitivity errorbars')

    figure(figsize=(12,8))
    gcf().add_axes([0.1,0.2,0.85,0.75])


    means=np.empty((len(outrange),2*ndim))
    stds=np.empty((len(outrange),2*ndim))
    for i in range(len(outrange)):
        iout=outrange[i]
        print('======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') =========')

        results_sens=pick.load(open('results_sens_'+str(iout+1)+'.pk', 'rb'))
        for j in range(2*ndim):
            mean=results_sens['pcmi'][0][j][0]
            varfrac0=results_sens['pcmi'][2][j][0]
            std=mean/sqrt(varfrac0)

            means[i,j]=mean
            stds[i,j]=std


    for maintot in ['main', 'total']:
        if (maintot=='main'):
            maintot_ind=range(0,ndim)
        elif (maintot=='total'):
            maintot_ind=range(ndim,2*ndim)

        figure(figsize=(12,8))
        gcf().add_axes([0.1,0.2,0.85,0.75])


        for i in range(len(outrange)):
            iout=outrange[i]
            errorbar(range(1,ndim+1),means[i,maintot_ind],yerr=[stds[i,maintot_ind],stds[i,maintot_ind]], fmt='-o',markersize=3, label=outnames[iout])

        #yscale('log')
        xlim([0.5,ndim+0.5])
        xticks(range(1,ndim+1),pnames,rotation=45)
        legend()
        savefig('senserb_par_'+maintot+'.eps')
        clf()

        figure(figsize=(12,8))
        gcf().add_axes([0.1,0.2,0.85,0.75])
        for j in maintot_ind:
            errorbar(outrange,means[:,j],yerr=[stds[:,j],stds[:,j]], fmt='-o',markersize=3, label=pnames[j%ndim])
        xlim([outrange[0]-0.5,outrange[-1]+0.5])
        xticks(outrange,outnames,rotation=45)
        legend()
        savefig('senserb_out_'+maintot+'.eps')
        clf()


else:
    print("plot_type not recognized. Exiting.")
    sys.exit()
