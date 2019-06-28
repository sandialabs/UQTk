#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version 3.0.4
#                     Copyright (2017) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#===================================================================================== 

import os
import shutil
import sys
try:
    import numpy as np
except ImportError:
    print "Numpy was not found. "

try:
    import matplotlib
except ImportError:
    print "Matplotlib was not found. "
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab


sys.path.append(os.environ['UQTK_SRC'])
import PyUQTk.plotting.surrogate as ss
import PyUQTk.plotting.inout as io

from pylab import *

import cPickle as pick

rc('legend',loc='best', fontsize=22)
rc('lines', linewidth=2, color='r')
rc('axes',linewidth=3,grid=True,labelsize=28)
rc('xtick',labelsize=20)
rc('ytick',labelsize=20)

    
#############################################################
#############################################################

plotid=(sys.argv[1])

results=pick.load(open('results.pk', 'rb'))
nout=results['training'][2].shape[1]
ndim=results['sens'][0].shape[1]
outrange=results['outs']



if os.path.exists('pnames.txt'):
    with open('pnames.txt') as f:
        pnames = f.read().splitlines()
        assert(len(pnames)==ndim)
else:
    pnames=['Param # '+str(i) for i in range(1,ndim+1)]

if os.path.exists('outnames.txt'):
    with open('outnames.txt') as f:
        outnames = f.read().splitlines()
        assert(len(outnames)==nout)
else:
    outnames=['Output # '+str(i) for i in range(1,nout+1)]


if(plotid=='sens'):
    allsens_main,allsens_total,allsens_joint=results['sens']
    
    sensmode=sys.argv[2]
    if sensmode=='main':
        sensdata=allsens_main
    elif sensmode=='total':
        sensdata=allsens_total
    
    pars=range(ndim) 
    cases=range(nout)
    
    np.savetxt('allsens_'+sensmode+'.dat',sensdata)


    print '==================================================='
    print "Plotting "+sensmode+" sensitivities across all output QoIs (bar plot)"  
    ss.plot_sens(sensdata[outrange],pars,cases,vis="bar",reverse=False,par_labels=pnames,case_labels=outnames,ncol=4,grid_show=False, xlbl='', showplot=False)
    os.system('mv sens.eps sens_'+sensmode+'.eps')
    print 'Created file sens_'+sensmode+'.eps'

elif(plotid=='senscirc'):
    allsens_main,allsens_total,allsens_joint=results['sens']
 
    
    npar=allsens_main[0].shape[0]
    allsens_main_ave=np.zeros((npar,))
    allsens_joint_ave=np.zeros((npar,npar))
    print 'Plotting main sand joint sensitivities (circular plots)'
    for i in range(len(outrange)):
        iout=outrange[i]
        print '======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') ========='
        ss.plot_senscirc(outnames[iout], allsens_main[i],allsens_joint[i],pnames,showplot=False)
        os.system('mv senscirc.eps senscirc_output_'+str(iout+1)+'.eps')
        print 'Created file senscirc_output_'+str(iout+1)+'.eps'

        allsens_main_ave+=allsens_main[i]/npar
        allsens_joint_ave+=allsens_joint[i]/npar

    print 'Plotting averaged main sand joint sensitivities (circular plots)'
    ss.plot_senscirc('', allsens_main_ave,allsens_joint_ave,pnames,showplot=False)
    os.system('mv senscirc.eps senscirc_ave.eps')
    print 'Created file senscirc_ave.eps'


elif(plotid=='sensmat'):
    allsens_main,allsens_total,allsens_joint=results['sens']

    sensmode=sys.argv[2]
    if sensmode=='main':
        sensdata=allsens_main
    elif sensmode=='total':
        sensdata=allsens_total
    
    
    pars=range(ndim) 
    cases=range(nout)

    print '==================================================='
    print "Plotting "+sensmode+" sensitivities across all output QoIs (matrix plot)"  
    ss.plot_sensmat(sensdata[outrange],pars,cases,vis="bar",reverse=False,par_labels=pnames,case_labels=outnames,showplot=False)
    os.system('mv sensmat.eps sensmat_'+sensmode+'.eps')
    print 'Created file sensmat_'+sensmode+'.eps'


elif(plotid=='dm'):
    
    # Parse the arguments
    trvals=sys.argv[2:] 

    for i in range(len(outrange)):
        iout=outrange[i]
        axes_labels=['Model ('+outnames[iout]+')','Polynomial Surrogate']
        print '======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') ========='
        datas=[]
        models=[]
        errorbars=[]
        labels=[]
        for trval in trvals:
            print "Plotting model-vs-surrogate at "+trval+" points"
            datas.append(results[trval][2][:,iout])
            models.append(results[trval][3][:,i])
            errorbars.append([results[trval][4][:,i],results[trval][4][:,i]])
            labels.append(trval+' points')
        
        ss.plot_dm(datas,models,errorbars,labels,axes_labels,showplot=False)
        os.system('mv dm.eps dm_output_'+str(iout+1)+'.eps')
        print 'Created file dm_output_'+str(iout+1)+'.eps'

elif(plotid=='idm'):
    
    # Parse the arguments
    trval=sys.argv[2] 

    for i in range(len(outrange)):
        iout=outrange[i]
        #axes_labels=['Model ('+outnames[iout]+')','Polynomial Surrogate']
        print '======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') ========='
        
        print "Plotting runId-vs-model and runId-vs-surrogate at "+trval+" points"
        data=results[trval][2][:,iout]
        model=results[trval][3][:,i]
        errbar=[results[trval][4][:,i],results[trval][4][:,i]]
        #label='trval+' points

        ss.plot_idm(data,model,errbar,sort='none',figname='idm.eps')

        os.system('mv idm.eps idm_output_'+trval+'_'+str(iout+1)+'.eps')
        print 'Created file idm_output_'+trval+'_'+str(iout+1)+'.eps'


elif(plotid=='mindex'):

    print 'Plotting multiindices'

    for i in range(len(outrange)):
        iout=outrange[i]
        print '======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') ========='

        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        varfrac=results['pcmi'][2][iout]

        np.savetxt("pcf_output_"+str(iout+1)+".dat",pcf)
        np.savetxt("mindex_output_"+str(iout+1)+".dat",mindex,fmt='%d')
        ss.plot_mindex(mindex[1:,:],varfrac[1:],str(outnames[iout]),showplot=False)
        os.system('mv mindex.eps mindex_output_'+str(iout+1)+'.eps')
        print 'Created file mindex_output_'+str(iout+1)+'.eps'

elif(plotid=='micf'):
    if (ndim>3):
        print "plot_micf utility only visualizes 2d or 3d multiindices. Exiting."
        sys.exit()

    for i in range(len(outrange)):
        iout=outrange[i]
        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        np.savetxt("pcf_output_"+str(iout+1)+".dat",pcf)
        np.savetxt("mindex_output_"+str(iout+1)+".dat",mindex,fmt='%d')
        ss.plot_micf(mindex,cfs=pcf,showplot=False)
        os.system('mv micf.eps micf_output_'+str(iout+1)+'.eps')
        print 'Created file micf_output_'+str(iout+1)+'.eps'

elif(plotid=='pdf'):
    nsam=100000

    print 'Computing and plotting output PDFs'

    for i in range(len(outrange)):
        iout=outrange[i]    
        print '======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') ========='

        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]
        pctype=results['pcmi'][-1]
        order= np.amax(mindex)

        custom_xlabel='Output QoI ('+outnames[iout]+')'

        ss.plot_pcpdf(pctype,mindex,pcf,nsam,custom_xlabel,showplot=False)
        os.system('mv pcdens.eps pcdens_output_'+str(iout+1)+'.eps')
        print 'Created file pcdens_output_'+str(iout+1)+'.eps'

elif(plotid=='senserb1'):
    sensord=1
    senssam=1000
    print 'Computing sensitivity errorbars'

    means=np.empty((len(outrange),ndim))
    stds=np.empty((len(outrange),ndim))
    for i in range(len(outrange)):
        iout=outrange[i]    
        print '======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') ========='

        pcf=results['pcmi'][0][iout]
        mindex=results['pcmi'][1][iout]

        cfcov=results['pcmi'][3][iout]
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


elif(plotid=='senserb2'):
    
    print 'Plotting sensitivity errorbars'

    figure(figsize=(12,8))
    gcf().add_axes([0.1,0.2,0.85,0.75])


    means=np.empty((len(outrange),2*ndim))
    stds=np.empty((len(outrange),2*ndim))
    for i in range(len(outrange)):
        iout=outrange[i]    
        print '======== Output # '+str(iout+1) +' (QoI ' +outnames[iout]+') ========='

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
    print "plotid not recognized. Exiting."
    sys.exit()
