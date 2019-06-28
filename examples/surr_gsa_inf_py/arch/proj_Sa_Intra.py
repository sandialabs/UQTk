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
from dataTypeConvert import *

#############################################################
#############################################################
class bscSA(object):
  

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


#############################################################
  def initPrms(self):
    self.dir = os.getcwd()
    self.dataDump = os.path.join(self.dir,"dataDumpProj")
    if not os.path.exists(self.dataDump):
      os.makedirs(self.dataDump)
    
    # Process parameters
    logp = np.log10(self.prmsVal)
    self.ptrain = np.array([[logp-0.25],[logp+0.25]])
    

#############################################################
  def genQuadPnts(self):
    # get quad points and weights
    x_uq = mkBlnkDbl2D();
    w_uq = mkBlnkDbl1D();
            
    print "Create an instance of Quad class"
    level = 1
    q = uqtkquad.Quad('LU','sparse',self.ndim,level,0,1)
    q.SetRule()
    q.GetRule(x_uq,w_uq)
    self.q = q;

    print "Number of quad points is " + str(len(x_uq))
    self.w_np = uqtk2NumpyDbl(w_uq)
    self.x_np = uqtk2NumpyDbl(x_uq)
    #self.x_np = self.x_np*self.prmsVal
    
    self.saveDat(self.x_np,'x.dat')
    self.saveDat(self.w_np,'w.dat')
  
    
############################################################# 
  def runModel(self):
    msk = np.arange(1,self.ndim+1)
    msk = 1.0/(msk);#*msk
   
    y   = np.dot(self.x_np,msk);
    self.y_np = y;
    self.saveDat(self.y_np,'y.dat')


############################################################# 
  def runPCE(self):
    print "Instantiate PCE object"
    self.nord = 1;
    pcmodel = uqtkpce.PCSet('NISPnoq',self.nord,self.ndim,'LEG')
    # set quad rule for pc model
    pcmodel.SetQuadRule(self.q)

    # Get the multiindex for postprocessing
    mIdx_uq  = mkBlnkInt2D();
    pcmodel.GetMultiIndex(mIdx_uq);
    self.mIdx_np = uqtk2NumpyDbl(mIdx_uq)
    self.saveDat(self.mIdx_np,'mIdex.dat')

    # get the coefficients using the quadrature rule to calculate the projections
    nPCT = pcmodel.GetNumberPCTerms();
    ck_uq = mkSzDbl1D(nPCT);
    y_uq = mkSetDbl1D(self.y_np);
    pcmodel.GalerkProjection(y_uq,ck_uq);
    self.ck_np = uqtk2NumpyDbl(ck_uq)
    self.saveDat(self.ck_np,'ck.dat')
   
    self.pcModel = pcmodel;
    self.mIdx_uq = mIdx_uq;
    self.ck_uq = ck_uq;


############################################################# 
  def runSens(self,nm):
    pcmodel = self.pcModel;
    
    # compute main sensitivities
    mainsens_uq = mkSzDbl1D(self.ndim)
    pcmodel.ComputeMainSens(self.ck_uq,mainsens_uq)

    # compute total sensitivity
    totalsens_uq = mkSzDbl1D(self.ndim)
    pcmodel.ComputeTotSens(self.ck_uq,totalsens_uq)

    #compute joint sensitivity
    jointsens_uq = mkSzDbl2D(self.ndim,self.ndim)
    pcmodel.ComputeJointSens(self.ck_uq,jointsens_uq)
  
    self.mainsens_np  = uqtk2NumpyDbl(mainsens_uq)
    self.totalsens_np = uqtk2NumpyDbl(totalsens_uq)
    self.jointsens_np = uqtk2NumpyDbl(jointsens_uq)
    self.saveDat(self.mainsens_np ,'sens'+ nm + 'Main.dat')
    self.saveDat(self.totalsens_np,'sens'+ nm + 'Total.dat')
    self.saveDat(self.jointsens_np,'sens'+ nm + 'Joint.dat')
   

#############################################################
  def main(self):
    #Processes Parameters
    self.nBCS = 2;
    self.rdPrmsFle();
    self.initPrms();
    print "Parameters Processsed \n"
    
    #Generate Quad Points
    self.genQuadPnts()
    print "Quad Points Generated \n"
    
    #run model with 
    self.runModel()
    print "Models Evalulated \n"
    
    #run polynomial chaos expansion
    self.runPCE()
    print "Polynomial Chaos Expansion Built \n"
    
    #run sensitivity analysis
    self.runSens('Std')
    print "Projection Sensitivity Analysis \n"
    

#############################################################
if __name__ == "__main__":
    bscSA().main()
    
    
    
    
    
    
    
    
    
    
    
    
