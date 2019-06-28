import os
import numpy as npy
import matplotlib.pyplot as plt

global xdata, ydata, sigma

def loglik(theta,xdata,ydata,sig2):
    nord = theta.shape[1]-1
    nspl = theta.shape[0]
    npts = xdata.shape[0]
    llik = npy.zeros(nspl)
    for j in range(nspl):
      for i in range(npts):
        if nord>0:
            peval = theta[j,nord]
            for k in range(nord-1,-1,-1):
                peval = peval*xdata[i]+theta[j,k]
        else:
            peval = theta[j,0]
        llik[j] -= ((peval-ydata[i])**2/(2*sig2)+0.5*npy.log(2.0*npy.pi*sig2))
    #llik[npy.where(llik<npy.log(1.e-80))]=npy.log(1.e-80)
    return llik


#------------------------------------------
#  load data
#------------------------------------------
npts=11
sig2=0.01
ydata=npy.genfromtxt("yd_n11_s0.01.dat")
xdata=npy.linspace(-1,1,npts)

#------------------------------------------
#  load theta
#------------------------------------------
fmcmc="mcmcstates_1.dat"
if os.stat(fmcmc).st_size == 0:
  quit()

theta=npy.genfromtxt(fmcmc)
if len(theta.shape)==1:
  theta = npy.array([theta])

llik=loglik(theta,xdata,ydata,sig2);
npy.savetxt('tmcmc_ll.dat', llik, fmt='%.18e', delimiter=' ', newline='\n')


