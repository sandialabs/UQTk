import PyUQTk.uqtkarray as uqtkarray
import PyUQTk.tools as uqtktools
import numpy as np

#~# Double Arrays
#############################################################
def mkSetDbl1D(v_np):
  v_uqtk = uqtkarray.dblArray1D(v_np.shape[0])
  v_uqtk.setnpdblArray(v_np)
  return v_uqtk


#############################################################
def mkSetDbl2D(v_np):
  v_uqtk = uqtkarray.dblArray2D(v_np.shape[0],v_np.shape[1])
  v_uqtk.setnpdblArray(v_np)
  return v_uqtk


#############################################################
def mkSzDbl1D(dim1):
  return uqtkarray.dblArray1D(dim1)


#############################################################
def mkSzDbl2D(dim1, dim2):
  return uqtkarray.dblArray2D(dim1, dim2)


#############################################################
def mkBlnkDbl1D():
  return uqtkarray.dblArray1D()


#############################################################
def mkBlnkDbl2D():
  return uqtkarray.dblArray2D()


#############################################################
def uqtk2NumpyDbl(v_uqtk):
  v_np = uqtkarray.uqtk2numpy(v_uqtk)
  return v_np



#~# Integer Arrays
#############################################################
def mkSetInt1D(v_np):
  v_uqtk = uqtkarray.intArray1D(v_np.shape[0])
  v_uqtk.setnpdblArray(v_np)
  return v_uqtk


#############################################################
def mkSetInt2D(v_np):
  v_uqtk = uqtkarray.intArray2D(v_np.shape[0],v_np.shape[1])
  v_uqtk.setnpdblArray(v_np)
  return v_uqtk


#############################################################
def mkSzInt1D(dim1):
  return uqtkarray.intArray1D(dim1)


#############################################################
def mkSzInt2D(dim1, dim2):
  return uqtkarray.intArray2D(dim1, dim2)


#############################################################
def mkBlnkInt1D():
  return uqtkarray.intArray1D()


#############################################################
def mkBlnkInt2D():
  return uqtkarray.intArray2D()




#~# Double Reference
#############################################################
def mkSetDblRef(v_np):
  v_uqtk  = uqtktools.new_doublep()
  uqtktools.doublep_assign(v_uqtk,v_np)
  return v_uqtk


#############################################################
def mkDblRef():
  v_uqtk  = uqtktools.new_doublep()
  return v_uqtk





