"Pybind File for the Array module creation"
from PYB11Generator import *

PYB11includes = ["arrayio.h","arraytools.h","stdlib.h","stdio.h","math.h","assert.h","sstream","fstream","iomanip","ftndefs.h","gen_defs.h","depblas.h","deplapack.h"]

PYB11namespaces = ["std"]

def read_datafile():
    "Read a datafile from filename in a matrix form and store it in a 2D array"
    return "void"

def read_datafileVS():
    "Read a datafile from filename in a vector form and store it in an array data of typename T"
    return "void"

def read_datafile_1d():
    "Read a datafile from filename in a vector form and store it in a 1d array data of typename T"
    return "void"

def write_datafile_size():
    "Write the contents of a 2d array data of typename T to file filename in a matrix form"
    return "void"

def write_datafile():
    "Write the contents of an array data of typename T to file filename in a matrix form"
    return "void"

def write_datafile_1d():
    "Write the contents of a 1d array data of typename T to file filename in a vector form"
    return "void"

def array1Dto2D():
    "Store a given 1d array in a 2d array with a single second dimension"
    return "void"

def array2Dto1D():
    "Store a given 2d array in a 1d array with a single second dimension"
    return "void"

def paste():
    "Paste two arrays of same size into a single array with second dimension equal to two"
    return "void"

def generate_multigrid():
    "Generates multigrid as a cartesian product of each column of grid"
    return "void"

def merge(x = "Array2D<double>",y = "Array2D<double>",xy = "Array2D<double>"):
    "Merges 2D double arrays"
    return "void"

@PYB11pycppname("merge")
def merge1(x = "Array1D<double>",y = "Array1D<double>",xy = "Array1D<double>"):
    "Merges 1D double arrays"
    return "void"

@PYB11pycppname("merge")
def merge2(x = "Array1D<int>",y = "Array1D<int>",xy = "Array1D<int>"):
    "Merges 1D int arrays"
    return "void"

def append(x = "Array1D<double>",y = "Array1D<double>"):
    "Append one array to another (double format)"
    return "void"

@PYB11pycppname("append")
def append1(x = "Array1D<int>",y = "Array1D<int>"):
    "Append one array to another (int format)"
    return "void"

def transpose():
    "Transpose a 2d double or int array x and return the result in xt"
    return "void"

def flatten():
    "Flatten a 2d array into a 1d array"
    return "void"

def fold_1dto2d_rowfirst():
    "Fold a 1d array into a 2d array (double format), row first"
    return "void"

def fold_1dto2d_colfirst():
    "Fold a 1d array into a 2d array (double format), column first"
    return "void"

def swap(arr = "Array1D<double>",i ="int", j = "int"):
    "Swap i-th and j-th elements of the array arr"
    return "void"

@PYB11pycppname("swap")
def swap1(arr = "Array2D<double>",i ="int", j = "int"):
    "Swap i-th and j-th row of the 2d array arr"
    return "void"

def access():
    "Access element of an array"
    return "void"

def getRow():
    "Retrieves row 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'"
    return "void"

def getCol():
    "Retrieves column 'k' from 2D array 'arr2d' and returns it in 1D array 'arr1d'"
    return "void"

def addVal():
    "Adds 'val' to all elements of 1D array arr1d (double or int)"
    return "void"

def subVector():
    "Extracts from 'vector', elements corresponding to indices 'ind' and returns them in 'subvector' (double or int)"
    return "void"

def subMatrix_row():
    "Extracts from 'matrix' rows corresponding to indices 'ind' and returns them in 'submatrix' (double or int)"
    return "void"

def subMatrix_col():
    "Extracts from 'matrix' columns corresponding to indices 'ind' and returns them in 'submatrix' (double or int)"
    return "void"

def matPvec():
    "Adds scaled row or column to all rows / columns of a matrix (double or int)"
    return "void"

def maxVal():
    "Returns maximum value in 'vector' and its location in *indx (double or int)"
    return "void"

def setDiff():
    "Returns in C elements of A that are not in B; C is sorted in ascending order"
    return "void"

def setdiff_s():
    "Returns in C elements of A that are not in B; C is sorted in ascending order. Assumes A is sorted"
    return "void"

def shell_sort(a = "int *",n = "int"):
    "Sorts integer array"
    return "void"

@PYB11pycppname("shell_sort")
def shell_sort1(array = "Array1D<int>"):
    "Sorts integer array in ascending order"
    return "void"

@PYB11pycppname("shell_sort")
def shell_sort2(array = "Array1D<double>"):
    "Sorts double array in ascending order"
    return "void"

def shell_sort_col():
    "Sorts double array in ascending order according to a given column"
    return "void"

def shell_sort_all():
    "Sorts double array in ascending order according to first column, then second column breaks the tie, and so on"
    return "void"

def quicksort3(arr = "Array1D<double>",l = "int", r = "int"):
    "Quick-sort with 3-way partitioning of array between indices l and r"
    return "void"

@PYB11pycppname("quicksort3")
def quicksort31(arr = "Array2D<double>",l = "int", r = "int"):
    "Quick-sort with 3-way partitioning of 2d array between indices l and r and sorting is done comparing rows (by first element, then by second, etc...)"
    return "void"

@PYB11pycppname("quicksort3")
def quicksort32(arr = "Array2D<double>",l = "int", r = "int",col = "int"):
    "Quick-sort with 3-way partitioning of 2d array between indices l and r, according to column col"
    return "void"

def intersect(A = "Array1D<int>", B = "Array1D<int>",C = "Array1D<int>", iA = "Array1D<int>",iB = "Array1D<int>"):
    "Finds common entries in 1D arrays 'A' and 'B' and returns them in 'C', sorted in ascending order. It also returns the original locations of these entries in 1D arrays 'iA' and 'iB', respectively"
    return "void"

@PYB11pycppname("intersect")
def intersect1(A = "Array1D<int>", B = "Array1D<int>",C = "Array1D<int>"):
    "Finds common entries in 1D arrays 'A' and 'B' and returns them in 'C', sorted in ascending order."
    return "void"

def find():
    "Return list of indices corresponding to elements of 1D array theta that are: larger ( type='gt' ), larger or equal ( type='ge' ), smaller ( type='lt' ), smaller or equal ( type='le' ) than lmbda"
    return "void"

def prodAlphaMatVec():
    "Implements y = a A x"
    return "void"

def prodAlphaMatTVec():
    "Implements y = a A^T x"
    return "void"

def prodAlphaMatMat():
    "Implements C = a A B"
    return "void"

def prodAlphaMatTMat():
    "Implements C = a A^T B"
    return "void"

def addVecAlphaVecPow():
    "Implements x = x + a y**ip"
    return "void"

def prod_vecTmatvec():
    "Returns a^T B c"
    return "double"

def MatTMat():
    "Returns A^T A"
    return "Array2D<double>"

def delRow():
    "Deletes a row from a matrix"
    return "void"

def delCol(A = "Array2D<T>",icol = "int"):
    "Deletes a column from a matrix"
    return "void"

@PYB11pycppname("delCol")
def delCol1(A = "Array1D<T>",icol = "int"):
    "Deletes an element from a matrix"
    return "void"

def paddMatRow():
    "Padds 2D array 'A' with the row 'x'"
    return "void"

def paddMatCol)():
    "Padds 2D array 'A' with the column 'x'"
    return "void"

def paddMatColScal():
    "Padds square matrix A with a row and column x, and adds an element A_{n+1,n+1} to obtain a larger square matrix"
    return "void"

def is_equal(a="Array1D<int>",b="Array1D<int>"):
    "Checks if two 1d int arrays are equal"
    return "bool"

@PYB11pycppname("is_equal")
def is_equal1(a="Array1D<double>",b="Array1D<double>"):
    "Checks if two 1d double arrays are equal"
    return "bool"

def is_less(a="Array1D<int>",b="Array1D<int>"):
    "Checks if one 1d int array is less than another (by first element, then by second, etc...)"
    return "bool"

@PYB11pycppname("is_less")
def is_less1(a="Array1D<double>",b="Array1D<double>"):
    "Checks if one 1d double array is less than another (by first element, then by second, etc...)"
    return "bool"

def vecIsInArray():
    "Checks if vec matches with any of the rows of a 2d array"
    return "int"

def select_kth():
    "Select the k-th smallest element of an array arr"
    return "double"

def logdeterm():
    "Log-determinant of a real symmetric positive-definite matrix"
    return "double"

def trace():
    "Computes trace of a matrix"
    return "double"

def evalLogMVN():
    "Evaluates the natural logarithm of a multivariate normal distribution"
    return "double"

def diag():
    "Returns a diagonal matrix with a given diagonal"
    return "Array2D<double>"

def copy():
    "Returns a copy of a 1D array"
    return "Array1D<double>"

def copy():
    "Returns a copy of a 2D array"
    return "Array2D<double>"

def mtxdel():
    "Deletes matrix columns or rows"
    return "Array2D<double>"

def add():
    "Adds two vectors"
    return "Array1D<double>"

def add():
    "Adds two matrices of the same size"
    return "Array2D<double>"

def addinplace(x="Array2D<double>",y="Array2D<double>"):
    "Add two matrices of the same size"
    return "void"

@PYB11pycppname("addinplace")
def addinplace1(x="Array1D<double>",y="Array1D<double>"):
    "Add two vectors of the same size"
    return "void"

def subtract():
    "Subtract two vectors"
    return "Array1D<double>"

def subtract():
    "Subtract two matrices of the same size"
    return "Array2D<double>"

def subtractinplace(x="Array2D<double>",y="Array2D<double>"):
    "Subtract two matrices of the same size"
    return "void"

@PYB11pycppname("subtractinplace")
def subtractinplace1(x="Array1D<double>",y="Array1D<double>"):
    "Subtract two vectors of the same size"
    return "void"

def scale():
    "Multiply 1D Array by a double"
    return "Array1D<double>"

def scale():
    "Multiply 2D Array by a double"
    return "Array2D<double>"

def scaleinplace(x="Array1D<double>",alpha="double"):
    "multiply Array1D by double, in place"
    return "void"

@PYB11pycppname("scaleinplace")
def scaleinplace1(x="Array2D<double>",alpha="double"):
    "multiply 2d Array by double, in place"
    return "void"

@PYB11pycppname("scaleinplace")
def scaleinplace2(x="Array2D<double>",alpha="int"):
    "multiply 2d Array by int, in place"
    return "void"

def dotmult():
    "Returns the elementwise multiplication of two 2D Arrays"
    return "Array2D<double>"

def dotmult():
    "Returns the elementwise multiplication of two 1D Arrays"
    return "Array1D<double>"

def dotdivide():
    "Returns the elementwise division of two 2D Arrays"
    return "Array2D<double>"

def dotdivide():
    "Returns the elementwise division of two 1D Arrays"
    return "Array1D<double>"

def norm():
    "Get norm of array 1d"
    return "double"

def dist_sq():
    "Returns weighted vector distance-squared"
    return "double"

def Trans():
    "Get transpose of 2d array"
    return "Array2D<double>"

def dot():
    "get dot product between two vectors"
    return "double"

def dot(A="Array2D<double>",x="Array1D<double>"):
    "get matrix vector product"
    return "Array1D<double>"

@PYB11pycppname("dot")
def dot1(x="Array1D<double>",A="Array2D<double>"):
    "get matrix vector product"
    return "Array1D<double>"

def dot():
    "get matrix matric product"
    return "Array2D<double>"

def dotT():
    "get matrix^T matrix product"
    return "Array2D<double>"

def INV():
    "Inverse of a real square matrix"
    return "Array2D<double>"

def AinvH():
    "Solving AX=H where A is real, symmetric and positive definite"
    return "Array2D<double>"

def Ainvb():
    "Solving Ax=b where A is real, symmetric and positive definite"
    return "Array1D<double>"

def LSTSQ():
    "Least squares solution for overdetermined system"
    return "void"

def QR():
    "QR factorization"
    return "void"

def SVD():
    "SVD calculation"
    return "void"

def printarray(x="Array1D<double>"):
    "Print array to screen"
    return "void"

@PYB11pycppname("printarray")
def printarray1(x="Array1D<int>"):
    "Print array to screen"
    return "void"

@PYB11pycppname("printarray")
def printarray2(x="Array2D<int>"):
    "Print array to screen"
    return "void"

@PYB11pycppname("printarray")
def printarray3(x="Array2D<double>"):
    "Print array to screen"
    return "void"
