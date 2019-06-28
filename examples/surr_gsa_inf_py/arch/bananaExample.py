import numpy as npy
import scipy.stats
import scipy.linalg
import math
import matplotlib.pyplot as plt
from PyUQTk.inference.mcmc import dram

#############################################################
#############################################################
class bananaExample(object):

#############################################################
  def postrrBanana(self,spl,postinfo):
    # Rosenbrock Equation 
    # -((A-spl[0])**2+B*(spl[1]-spl[0]**2)**2)
    # a = 1
    # b = 100
    return -((1-spl[0])**2+100*(spl[1]-spl[0]**2)**2)


#############################################################
  def main(self):
    # defind variables
    cini  = npy.array([-1.0,- 4.0])
    spllo = npy.array([-4.0,-12.0])
    splhi = npy.array([ 4.0,  2.0])
    cvini = npy.array([[0.1,0.0],[0.0,0.1]])
    opts = {'method':'am','nsteps':10000,'nburn':1000,'nadapt':100,'nfinal':10000000,'inicov':cvini,'coveps':1.e-10,'burnsc':5,'ndr':3,'drscale':[5,4,3],'spllo':spllo,'splhi':splhi}
    lpinfo = {'afac':[1.0,1.0],'cov': npy.array([[1,0.9],[0.9,1]]),'mu':npy.array([0.0,0.0])}

    #Running MCMC tests
    print "Running amcmc test"
    sol = dram(opts,cini,self.postrrBanana,lpinfo)
    
    #plot results
    plt.plot(sol[0][1000:,0],sol[0][1000:,1],'o');
    plt.ion()
    plt.show()
    raw_input("Press Enter to close")
    plt.close()


#############################################################
if __name__ == "__main__":
    bananaExample().main()
