import os
import string 
import numpy as np
import PyUQTk.uqtkarray
import PyUQTk.quad as uqtkquad
import PyUQTk.pce as uqtkpce
import PyUQTk.bcs as uqtkbcs
import PyUQTk.utils.multiindex as uqtkmlti
import PyUQTk.tools as uqtktools
#import PyUQTk.sens.gsalib as gsalib
from dataTypeConvert import *

#############################################################
#############################################################
class bscSA(object):
  

#############################################################
  def saveDat(self,v,nm):
    np.savetxt(os.path.join(self.dataDump,nm), v)
  

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
    x_uq    = mkSetDbl2D(self.x_np);
    y_uq    = mkSetDbl1D(self.y_np);
    eta     = 1.e-4;
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
      
      if (1.0*len(usd_np) / self.mIdx_np.shape[0]) > 0.5:
        print 'No reduction, skipping BCS'
        self.skp = self.skp + 1;
        print str(self.skp) + ' of ' + str(self.cnt) + 'Skipped so far'
        break;

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
  def runSens(self):
    pcmodel = self.pcModel;
    w_uq = self.wghts_uq;
    
    # compute main sensitivities
    mainsens_uq = mkSzDbl1D(self.ndim)
    pcmodel.ComputeMainSens(w_uq,mainsens_uq)
    self.mainsens_np  = uqtk2NumpyDbl(mainsens_uq)
    
    #compute joint sensitivity
    # jointsens_uq = mkSzDbl2D(self.ndim,self.ndim)
    # pcmodel.ComputeJointSens(w_uq,jointsens_uq)
    # self.jointsens_np = uqtk2NumpyDbl(jointsens_uq)

    # compute total sensitivity
    totalsens_uq = mkSzDbl1D(self.ndim)
    pcmodel.ComputeTotSens(w_uq,totalsens_uq)
    self.totalsens_np = uqtk2NumpyDbl(totalsens_uq)
    

#############################################################
  def main(self):
    print "\n\n Initalizing Pipeline \n" 
    self.nBCS = 3;
    self.nord = 1; 
    self.ndim = 254;
    self.nSpTm = 780;
    self.dir = os.getcwd()
    self.dataDump = os.path.join(self.dir,"dataDumpBCS")
    if not os.path.exists(self.dataDump):
      os.makedirs(self.dataDump)
    
    print 'Import xData and yData'
    dataDumpPrms = os.path.join(self.dir,"dataDumpPrms")
    self.x_np = np.genfromtxt(os.path.join(dataDumpPrms,"xData.csv"), delimiter=',')
    yFull = np.genfromtxt('yDataFull.csv', delimiter=',')
    print self.x_np.shape
    print yFull.shape
    
    #define sens matrix
    main  = np.zeros([self.ndim,self.nSpTm]);
    # joint = np.zeros([self.ndim,self.nSpTm]);
    total = np.zeros([self.ndim,self.nSpTm]);
    
    self.cnt = 1;
    self.skp = 0;
    for spTm in range(yFull.shape[1]):
      print 'running species per time point:' + str(self.cnt)
      self.y_np = yFull[:,spTm];

      #run BCS SA per species and time point 
      self.runPCE()
      self.runBCS()
      self.runSens()
      
      main[:,spTm]  = self.mainsens_np
      # joint[:,spTm] = self.jointsens_np
      total[:,spTm] = self.totalsens_np
      self.cnt = self.cnt + 1;
    
    self.saveDat(main ,'sensBcsMain.dat' )
    self.saveDat(total,'sensBcsTotal.dat')
    # self.saveDat(joint,'sensBcsJoint.dat')
    
    
#############################################################
if __name__ == "__main__":
    bscSA().main()
    
    
    
    
    
    
    
    
    
    
    
    
