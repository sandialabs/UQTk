import os
import string
import numpy as np
import pdb
import PyUQTk.uqtkarray
import PyUQTk.quad as uqtkquad
import PyUQTk.pce as uqtkpce
import PyUQTk.bcs as uqtkbcs
import PyUQTk.utils.multiindex as uqtkmlti
import PyUQTk.tools as uqtktools
import PyUQTk.sens.gsalib as gsalib
from DataTypeConvert import *

#############################################################
#############################################################
class BCSSA(object):


#############################################################
  def saveDat(self,v,nm):
    np.savetxt(os.path.join(self.dataDump,nm), v)


#############################################################
  def rdPrmsFle(self):
    #read parameters file
    f = open("params_hand.txt")
    lines = f.readlines()
    sz = len(lines);
    print  "\n" + "Parameters Read From File: " + str(sz)
    paramNM = []
    paramVL = []
    prmsVal = np.zeros([1,sz])

    for i in xrange(0,sz):
      tmpStr = lines[i];
      tmp = tmpStr.split('=', 1 )
      tmpNM,tmpVL = tmp[0],tmp[1]
      paramNM.append(tmpNM.strip(' \t\n\r,;'))
      paramVL.append(tmpVL.strip(' \t\n\r,;'))
      prmsVal[0,i] = float(paramVL[i])

    self.paramNM = paramNM;
    self.prmsVal = prmsVal;

    # Only Considering first 10 variables
    self.prmsVal = prmsVal[0,0:10];
    self.ndim = self.prmsVal.size
    print 'Only running SA on first ' + str(self.ndim) + ' vars'
    print self.prmsVal

#############################################################
  def initPrms(self):
    self.dir = os.getcwd()
    self.dataDump = os.path.join(self.dir,"dataDumpBCS")
    if not os.path.exists(self.dataDump):
      os.makedirs(self.dataDump)

    # Set range over which parameters are being sampled
    self.parRanges = np.array([(1.0 - self.halfwidth)*self.prmsVal,(1.0 + self.halfwidth)*self.prmsVal])


#############################################################
  def genLhsPnts(self):

    # lhsmdu module
    # lhs = lhsmdu.sample(self.ndim, self.nSpls)
    # self.x_np = np.dot(lhs.T,np.diag(self.prmsVal))
    # self.x_np += self.prmsVal

    # uqtk function
    lhs_uq = mkSzDbl2D(self.ndim, self.nSpls)
    uqtktools.generate_uniform_lhs(lhs_uq,0) # generates LHS sample points between 0 and 1


    print "Number of sample points is ", self.nSpls
    self.x_np = uqtk2NumpyDbl(lhs_uq)

    print self.x_np

    # Map LHS points from the standard [0,1] to the parameter values
    self.lam_np = self.parRanges[0,] + (self.parRanges[1,]-self.parRanges[0,])*self.x_np.T


    print self.lam_np

    self.saveDat(self.lam_np,'x.dat') # locations at which to evaluate forward model



#############################################################
  def runModel(self):
    """Forward model: \sum_{i=1}^{ndim} \lambda_i/i """
    msk = np.arange(1,self.ndim+1)
    msk = 1.0/(msk);#*msk

    y   = np.dot(self.lam_np,msk);
    self.y_np = np.array(y);
    self.saveDat(self.y_np,'y.dat')


#############################################################
  def runPCE(self):
    print "Instantiate PCE object"

    pcmodel = uqtkpce.PCSet('NISPnoq',self.nord,self.ndim,'LEG')
    mIdx_uq  = mkBlnkInt2D();
    pcmodel.GetMultiIndex(mIdx_uq);

    self.mIdx_np = uqtk2NumpyDbl(mIdx_uq)
    self.saveDat(self.mIdx_np,'mIdex.dat')

    self.pcModel = pcmodel;


#############################################################
  def runBCS(self):
    pcmodel = self.pcModel;
    x_uq    = mkSetDbl2D(self.x_np); # Should this be the original LHS points or the actual physical parameters that we fed into the model?
    y_uq    = mkSetDbl1D(self.y_np);
    eta     = 1.e-9;
    adpt    = 0;
    optml   = 1;
    scl     = 0.1;
    verb    = -1;

    for iters in range(1, self.nBCS+1):
      print 'Starting BCS ' + str(iters)
      self.nPCT   = pcmodel.GetNumberPCTerms();
      phi_uq      = mkBlnkDbl2D();
      midx_uq     = mkBlnkInt2D();
      sgm2_uq     = mkSetDblRef((np.std(self.y_np)**2)/1.0e6);
      lmbdInit_uq = mkBlnkDbl1D();
      wghts_uq    = mkBlnkDbl1D();
      usd_uq      = mkBlnkInt1D();
      errBrs_uq   = mkBlnkDbl1D();
      basis_uq    = mkBlnkDbl1D();
      alpha_uq    = mkBlnkDbl1D();
      lmbd_uq     = mkDblRef();
      newMIdx_uq  = mkBlnkInt2D();

      # evalulate pc model to get phi
      pcmodel.EvalBasisAtCustPts(x_uq, phi_uq)
      self.phi_np = uqtk2NumpyDbl(phi_uq)
      self.saveDat(self.phi_np,'phi_' + str(iters) + '.dat')

      # Perform BCS
      uqtkbcs.FastLaplace(phi_uq, y_uq, sgm2_uq, eta, lmbdInit_uq, adpt, optml, scl, verb, wghts_uq, usd_uq, errBrs_uq, basis_uq, alpha_uq, lmbd_uq)
      usd_np = uqtkarray.uqtk2numpy(usd_uq)
      print "BCS selected " + str(len(usd_np)) + " of " + str(self.mIdx_np.shape[0]) + " Multi Index"
      self.mIdx_np = self.mIdx_np[usd_np,:]

      # skip add front for last iteration of BCS
      if iters != self.nBCS:
        # Update Multi Index with BCS Results
        [mIdxNew_np, mIdxAdd_np, mIdxFrnt_np] = uqtkmlti.mi_addfront(self.mIdx_np)
        self.saveDat(mIdxNew_np,'nIdexNew_' + str(iters) + '.dat')

        # Redefine PCE with update Multi Index
        mIdxNew_uq = mkSetInt2D(mIdxNew_np);
        pcmodel = uqtkpce.PCSet('NISPnoq',mIdxNew_uq,'LEG')
        self.mIdx_np = uqtkarray.uqtk2numpy(mIdxNew_uq)
        print ''

    print ''
    mIdxFin_uq = mkSetInt2D(self.mIdx_np);
    pcmodel = uqtkpce.PCSet('NISPnoq',mIdxFin_uq,'LEG')
    print "Final mIdx size: " + str(self.mIdx_np.shape[0])
    self.pcModel = pcmodel;
    self.wghts_uq = wghts_uq;


#############################################################
  def runSens(self,nm):
    pcmodel = self.pcModel;
    w_uq = self.wghts_uq;
    # compute main sensitivities
    mainsens_uq = mkSzDbl1D(self.ndim)
    pcmodel.ComputeMainSens(w_uq,mainsens_uq)
    mainsens_np  = uqtk2NumpyDbl(mainsens_uq)
    self.saveDat(mainsens_np ,'sens'+ nm + 'Main.dat')

    # compute total sensitivity
    totalsens_uq = mkSzDbl1D(self.ndim)
    pcmodel.ComputeTotSens(w_uq,totalsens_uq)
    totalsens_np = uqtk2NumpyDbl(totalsens_uq)
    self.saveDat(totalsens_np,'sens'+ nm + 'Total.dat')

    #compute joint sensitivity
    jointsens_uq = mkSzDbl2D(self.ndim,self.ndim)
    pcmodel.ComputeJointSens(w_uq,jointsens_uq)
    jointsens_np = uqtk2NumpyDbl(jointsens_uq)
    self.saveDat(jointsens_np,'sens'+ nm + 'Joint.dat')


#############################################################
  def main(self):
    #Processes Parameters
    self.nSpls = 10; # Number of LHS samples
    self.nBCS = 3;   # Number of BCS iterations
    self.nord = 1;   # First order PCE since model is linear in parameters
    self.halfwidth = 0.2 # half the width of range around nominal parameter value that response surface
                         # will be created over expressed as ratio: i.e. if halfwidth is 0.2: sample in [0.8*a, 1.2*a]
    self.rdPrmsFle();
    self.initPrms();
    print "Parameters Processsed \n"

    #Generate Quad Points
    self.genLhsPnts()
    print "Latin Hypercube Points Generated \n"

    #run model with
    self.runModel()
    print "Models Evalulated \n"

    #run polynomial chaos expansion
    self.runPCE()
    print "Polynomial Chaos Expansion Built \n"

    #run basysian compress sensing
    self.runBCS()
    print "\n" + "Basysian Compress Sensing Completed \n"

    #run sensitivity analysis
    self.runSens('Bcs')
    print "BCS Sensitivity Analysis \n"


#############################################################
if __name__ == "__main__":
    BCSSA().main()
