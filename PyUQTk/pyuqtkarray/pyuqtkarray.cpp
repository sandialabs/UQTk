#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "Array1D.h"
#include "Array2D.h"
#include "arrayio.h"
#include "arraytools.h"
#include "ftndefs.h"
#include "gen_defs.h"
#include "depblas.h"
#include "deplapack.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include <sstream>
#include <fstream>
#include <iomanip>

namespace py=pybind11;

template <typename... Args>
using py_overload_cast = py::detail::overload_cast_impl<Args...>;

// Dummy function for read_datafileVS(Array2D<T> &data, const char *filename)
template <typename T>
void foo(Array2D<T> &data, const char *filename)
{

  ifstream in(filename);

  if(!in){
    printf("read_datafileVS() : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }

  string theLine="";

  // figure out number of lines and columns
  int nx, ny, ix = 0 ;
  while(in.good()){
    getline(in,theLine);

    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int    iy = 0 ;
    T      tmp    ;
    while( s >> tmp ) iy++ ;

    if ( ( ix > 0 ) && ( iy != ny ) )
    {
      printf("read_datafileVS() : Error at line %d !!!\n",ix+1) ;
      printf("                    no. of columns should be %d instead of %d\n",ny,iy) ;
      exit(1) ;
    }

    ny = iy ;

    ix++ ;

  }

  nx = ix ;

#ifdef VERBOSE
  printf("File \"%s\" contains %d rows and %d columns \n",filename,nx,ny) ;
#endif
  // Resize, goto beginning, and read again

  if ( ( (int) data.XSize() != nx ) || ( (int) data.YSize() != ny ) )
    data.Resize(nx,ny) ;

  //in.close() ;
  //in.open(filename);
  in.clear() ;
  in.seekg(0, ios::beg ) ;
  ix = 0 ;
  while( in.good() ){

    getline(in,theLine);

    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int    iy = 0 ;
    T      tmp ;
    while( s >> tmp ) {
      data(ix,iy)=tmp;
      iy++;
    }
    if ( iy != ny ) {
      printf("read_datafileVS() : Error in file \"%s\" \n",filename);
      printf("                    -> at line %d while reading %s; number of columns should be %d\n",
              ix+1, filename, ny);
      exit(1) ;
    }
    ix++;
  }
  if ( ix != nx ) {
    printf("read_datafileVS() : Error while reading \"%s\" -> number of rows should be %d\n", filename,nx) ;
    exit(1) ;
  }

  return ;

}
template void foo(Array2D<double> &data, const char *filename);
template void foo(Array2D<int>    &data, const char *filename);

// Dummy function for read_datafileVS(std::vector<T> &data, int &nrows, int &ncols, const char *filename)
template <typename T>
void fun(std::vector<T> &data, int &nrows, int &ncols, const char *filename)
{

  ifstream in(filename);

  if(!in){
    printf("read_datafileVS() : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }

  string theLine="";

  // figure out number of lines and columns
  int ix = 0 ;
  while(in.good()){
    getline(in,theLine);

    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int    iy = 0 ;
    T      tmp    ;
    while( s >> tmp ) iy++ ;

    if ( ( ix > 0 ) && ( iy != ncols ) )
    {
      printf("read_datafileVS() : Error at line %d !!!\n",ix+1) ;
      printf("                    no. of columns should be %d instead of %d\n",ncols,iy) ;
      exit(1) ;
    }

    ncols = iy ;

    ix++ ;

  }

  nrows = ix ;

#ifdef VERBOSE
  printf("File \"%s\" contains %d rows and %d columns \n",filename,nrows,ncols) ;
#endif
  // Resize, goto beginning, and read again

  data.resize(nrows*ncols,0.0);

  //in.close() ;
  //in.open(filename);
  in.clear() ;
  in.seekg(0, ios::beg ) ;
  ix = 0 ;
  while( in.good() ){

    getline(in,theLine);

    if ( theLine == "" ) break;
    if ( theLine.compare(0,1,"#") == 0 ) continue ;

    istringstream s(theLine);
    int    iy = 0 ;
    T      tmp ;
    while( s >> tmp ) {
      data[iy*nrows+ix]=tmp;
      iy++;
    }
    if ( iy != ncols ) {
      printf("read_datafileVS() : Error in file \"%s\" \n",filename);
      printf("                    -> at line %d while reading %s; number of columns should be %d\n",
              ix+1, filename, ncols);
      exit(1) ;
    }
    ix++;
  }
  if ( ix != nrows ) {
    printf("read_datafileVS() : Error while reading \"%s\" -> number of rows should be %d\n", filename,nrows) ;
    exit(1) ;
  }

  return ;

}
template void fun(std::vector<double> &data, int &nrows, int &ncols, const char *filename);
template void fun(std::vector<int>    &data, int &nrows, int &ncols, const char *filename);

// Dummy function for read_datafileVS(Array1D<T> &data, const char *filename)
template <typename T>
void dummy(Array1D<T> &data, const char *filename)
{

  data.Clear();
  ifstream in(filename);

  if(!in){
    printf("read_datafileVS() : the requested file %s does not exist\n",filename) ;
    exit(1) ;
  }

  string theLine="";
  int ix=0;

  while(in.good()){
    getline(in,theLine);

    if (theLine=="") break;

    istringstream s(theLine);
    T tmp;
    s >> tmp;
    data.PushBack(tmp);
    if (s>>tmp) {
      printf("Error at line %d while reading %s; number of columns should be 1\n", ix+1, filename);
      exit(1);
    }
    ix++;
  }

  return ;

}
template void dummy(Array1D<double> &data, const char *filename);
template void dummy(Array1D<int>    &data, const char *filename);

// Dummy function for write_datafile(const Array2D<T> &data, const char *filename)
template <typename T>
void rush(const Array2D<T> &data, const char *filename)
{

  int nx=data.XSize();
  int ny=data.YSize();

  FILE* f_out;
  if(!(f_out = fopen(filename,"w"))){
    printf("write_datafile: could not open file '%s'\n",filename);
    exit(1);
  }

  if ( typeid(T) == typeid(int) )
    for(int ix = 0 ; ix < nx ; ix++){
      for(int iy = 0 ; iy < ny ; iy++){
        fprintf(f_out, "%d ", data(ix,iy));
      }
      fprintf(f_out, "\n");
    }
  else if ( typeid(T) == typeid(double) )
    for(int ix = 0 ; ix < nx ; ix++){
      for(int iy = 0 ; iy < ny ; iy++){
        fprintf(f_out, "%24.16lg ", data(ix,iy));
      }
      fprintf(f_out, "\n");
    }
  else {
    printf("write_datafile: template not implemented\n");
    exit(1);
  }

  if(fclose(f_out)){
    printf("write_datafile: could not close file '%s'\n",filename);
    exit(1);
  }

#ifdef VERBOSE
  printf("Data written to '%s' in a matrix form [%d X %d]\n", filename,nx,ny);
#endif

 return ;

}
template void rush(const Array2D<double> &data, const char *filename);
template void rush(const Array2D<int>    &data, const char *filename);

// Dummy function for write_datafile(const Array2D<T> &data, const char *filename, const char *action)
template <typename T>
void stylebender(const Array2D<T> &data, const char *filename, const char *action)
{

  int nx=data.XSize();
  int ny=data.YSize();

  if ( ( string(action) != string("w") ) && ( string(action) != string("a") ) ) {
   printf("write_datafile: unknown file action '%s'\n",action);
   exit(1);
  }

  FILE* f_out;
  if(!(f_out = fopen(filename,action))){
    printf("write_datafile: could not open file '%s'\n",filename);
    exit(1);
  }

  if ( typeid(T) == typeid(int) )
    for(int ix = 0 ; ix < nx ; ix++){
      for(int iy = 0 ; iy < ny ; iy++){
        fprintf(f_out, "%d ", data(ix,iy));
      }
      fprintf(f_out, "\n");
    }
  else if ( typeid(T) == typeid(double) )
    for(int ix = 0 ; ix < nx ; ix++){
      for(int iy = 0 ; iy < ny ; iy++){
        fprintf(f_out, "%24.16lg ", data(ix,iy));
      }
      fprintf(f_out, "\n");
    }
  else {
    printf("write_datafile: template not implemented\n");
    exit(1);
  }

  if(fclose(f_out)){
    printf("write_datafile: could not close file '%s'\n",filename);
    exit(1);
  }

#ifdef VERBOSE
  printf("Data written to '%s' in a matrix form [%d X %d]\n", filename,nx,ny);
#endif

  return ;

}
template void stylebender(const Array2D<double> &data, const char *filename, const char *action);
template void stylebender(const Array2D<int>    &data, const char *filename, const char *action);

// Dummy function for write_datafile(const std::vector<T> &data, const int &nrows, const int &ncols, const char *storage, const char *filename, const char *action)
template <typename T>
void rowdy(const std::vector<T> &data, const int &nrows, const int &ncols, const char *storage, const char *filename, const char *action)
{

  if ( ( string(storage) != string("C") ) && ( string(storage) != string("R") ) ) {
   printf("write_datafile: unknown storage type '%s'\n",action);
   exit(1);
  }

  if ( ( string(action) != string("w") ) && ( string(action) != string("a") ) ) {
   printf("write_datafile: unknown file action '%s'\n",action);
   exit(1);
  }

  FILE *f_out;
  if(!(f_out = fopen(filename,action))){
    printf("write_datafile: could not open file '%s'\n",filename);
    exit(1);
  }

  if ( typeid(T) == typeid(int) ) {
    if ( string(storage) == string("C") ) {
      for(int ix = 0 ; ix < nrows ; ix++) {
        for(int iy = 0 ; iy < ncols ; iy++) {
          fprintf(f_out, "%d ", data[iy*nrows+ix]);
        }
        fprintf(f_out, "\n");
      }
    } else {
      for(int ix = 0 ; ix < nrows ; ix++) {
        for(int iy = 0 ; iy < ncols ; iy++) {
          fprintf(f_out, "%d ", data[ix*ncols+iy]);
        }
        fprintf(f_out, "\n");
      }
    }
  } // end of typeid int
  else if ( typeid(T) == typeid(double) ) {
    if ( string(storage) == string("C") ) {
      for(int ix = 0 ; ix < nrows ; ix++) {
        for(int iy = 0 ; iy < ncols ; iy++) {
          fprintf(f_out, "%24.16lg ", data[iy*nrows+ix]);
        }
        fprintf(f_out, "\n");
      }
    } else {
      for(int ix = 0 ; ix < nrows ; ix++) {
        for(int iy = 0 ; iy < ncols ; iy++) {
          fprintf(f_out, "%24.16lg ", data[ix*ncols+iy]);
        }
        fprintf(f_out, "\n");
      }
    }
  } // end of typeid double
  else {
    printf("write_datafile: template not implemented\n");
    exit(1);
  }

  if(fclose(f_out)){
    printf("write_datafile: could not close file '%s'\n",filename);
    exit(1);
  }

#ifdef VERBOSE
  printf("Data written to '%s' in a matrix form [%d X %d]\n", filename,nrows,ncols);
#endif

  return ;

}
template void rowdy(const std::vector<double> &data, const int &nrows, const int &ncols, const char *storage, const char *filename, const char *action);
template void rowdy(const std::vector<int>    &data, const int &nrows, const int &ncols, const char *storage, const char *filename, const char *action);

// Dummy methods for the dot function to work correctly
Array1D<double> method(Array2D<double>& A, Array1D<double>& x){
  int n=A.XSize();
  int m=A.YSize();

  // Size check
  CHECKEQ(m, (int) x.XSize());


  Array1D<double> y(n,0e0);

  char trans='n';
  double beta=0.e0;
  int xinc=1;
  int yinc=1;
  double alpha = 1;
  FTN_NAME(dgemv)(&trans, &n, &m, &alpha, A.GetArrayPointer(), &n, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );

  return y;
}
Array1D<double> swag(Array1D<double>& x, Array2D<double>& A){
  int n=A.XSize();
  int m=A.YSize();

  // Size check
  CHECKEQ(n, (int) x.XSize());


  Array1D<double> y(n,0e0);

  char trans='n';
  double beta=0.e0;
  int xinc=1;
  int yinc=1;
  double alpha = 1;
  FTN_NAME(dgemv)(&trans, &n, &m, &alpha, A.GetArrayPointer(), &n, x.GetArrayPointer(), &xinc,  &beta, y.GetArrayPointer(), &yinc );

  return y;
}

/*py::array_t<double> getnpdblArray(Array2D<double>& x){
  py::array_t<double> r_val(x.data_.size());
  auto buf = r_val.request();

  x.getnpdblArray((double * )buf.ptr);

  r_val.resize({x.XSize(),x.YSize()});

  return r_val;
}

py::array_t<int> getnpintArray(Array2D<int>& x){
  py::array_t<int> r_val(x.data_.size());
  auto buf = r_val.request();

  x.getnpintArray((long *)buf.ptr);

  r_val.resize({x.XSize(),x.YSize()});

  return r_val;
}*/

vector<int> getnpintArray(Array2D<int>& x){
  vector<int> a;
  int m = x.XSize();
  int n = x.YSize();
  for(int ip = 0;ip < m; ++ip){
    for(int idim = 0; idim < n; ++idim){
      a.push_back(x(ip,idim));
    }
  }
  return a;
}

vector<double> getnpdblArray(Array2D<double>& x){
  vector<double> a;
  int m = x.XSize();
  int n = x.YSize();
  for(int ip = 0;ip < m; ++ip){
    for(int idim = 0; idim < n; ++idim){
      a.push_back(x(ip,idim));
    }
  }
  return a;
}

py::array_t<int> getnpintArray(Array1D<int>& x){
  py::array_t<int> r_val(x.data_.size());
  auto buf = r_val.request();

  x.getnpintArray((long *)buf.ptr);

  return r_val;
}

py::array_t<double> getnpdblArray(Array1D<double>& x){
  py::array_t<double> r_val(x.data_.size());
  auto buf = r_val.request();

  x.getnpdblArray((double *)buf.ptr);

  return r_val;
}

void setnpdblArray(Array2D<double>& x,py::array_t<double> &inarray){
  /*  read input arrays buffer_info */
	auto buf1 = inarray.request();

  /*  variables */
	double *ptr1 = (double *) buf1.ptr;
	int n1 = buf1.shape[0];
  int n2 = buf1.shape[1];

	// Calling the class function
	x.setnpdblArray(ptr1,n1,n2);
}

void setnpintArray(Array2D<int>& x,py::array_t<int> &inarray){
  /*  read input arrays buffer_info */
	auto buf1 = inarray.request();

  /*  variables */
	long *ptr1 = (long *) buf1.ptr;
	int n1 = buf1.shape[0];
  int n2 = buf1.shape[1];

	// Calling the class function
	x.setnpintArray(ptr1,n1,n2);
}

PYBIND11_MODULE(_uqtkarray, m) {
    py::class_<Array1D<int>>(m, "intArray1D")
      .def(py::init<>())
      .def(py::init<const int&>())
      .def(py::init<const int&,const int&>())
      //.def("Assign", &Array1D<int>::operator=)
      .def(py::init<const Array1D<int> &>())
      .def("Clear",&Array1D<int>::Clear)
      .def("XSize",&Array1D<int>::XSize)
      .def("Length",&Array1D<int>::Length)
      .def("Resize",static_cast<void (Array1D<int>::*)(const int&)>(&Array1D<int>::Resize))
      .def("Resize",static_cast<void (Array1D<int>::*)(const int&,const int&)>(&Array1D<int>::Resize))
      .def("SetValue",&Array1D<int>::SetValue)
      .def("PushBack",&Array1D<int>::PushBack)
      .def("GetArrayPointer",&Array1D<int>::GetArrayPointer)
      .def("GetConstArrayPointer",&Array1D<int>::GetConstArrayPointer)
      .def("__getitem__", [](Array1D<int> a, int b) {
        return a(b);
      }, py::is_operator())
      .def("insert",static_cast<void (Array1D<int>::*)(Array1D<int>&,int)>(&Array1D<int>::insert))
      .def("insert",static_cast<void (Array1D<int>::*)(const int&,int)>(&Array1D<int>::insert))
      .def("erase",&Array1D<int>::erase)
      .def("DumpBinary",static_cast<void (Array1D<int>::*)(FILE*) const>(&Array1D<int>::DumpBinary))
      .def("DumpBinary",static_cast<void (Array1D<int>::*)(char*)>(&Array1D<int>::DumpBinary))
      .def("ReadBinary",static_cast<void (Array1D<int>::*)(FILE*)>(&Array1D<int>::ReadBinary))
      .def("ReadBinary",static_cast<void (Array1D<int>::*)(char*)>(&Array1D<int>::ReadBinary))
      .def("pyElement",&Array1D<int>::operator[])
      .def("ReadBinary4py",&Array1D<int>::ReadBinary4py)
      .def("DumpBinary4py",&Array1D<int>::DumpBinary4py)
      .def("setArray",&Array1D<int>::setArray)
      .def("setnpintArray",py::vectorize(&Array1D<int>::setnpintArray))
      //def("getnpintArray",py::vectorize(&Array1D<int>::getnpintArray))
      .def("flatten",&Array1D<int>::flatten)
      .def("type",&Array1D<int>::type)
      .def("shape",&Array1D<int>::shape)
      .def("assign",&Array1D<int>::assign)
      ;

      py::class_<Array1D<double>>(m,"dblArray1D")
        .def(py::init<>())
        .def(py::init<const int&>())
        .def(py::init<const int&,const double&>())
        //.def("Assign", &Array1D<double>::operator=)
        .def(py::init<const Array1D<double> &>())
        .def("Clear",&Array1D<double>::Clear)
        .def("XSize",&Array1D<double>::XSize)
        .def("Length",&Array1D<double>::Length)
        .def("Resize",static_cast<void (Array1D<double>::*)(const int&)>(&Array1D<double>::Resize))
        .def("Resize",static_cast<void (Array1D<double>::*)(const int&,const double&)>(&Array1D<double>::Resize))
        .def("SetValue",&Array1D<double>::SetValue)
        .def("PushBack",&Array1D<double>::PushBack)
        .def("GetArrayPointer",&Array1D<double>::GetArrayPointer)
        .def("GetConstArrayPointer",&Array1D<double>::GetConstArrayPointer)
        .def("__getitem__", [](Array1D<double> a, int b) {
          return a(b);
        }, py::is_operator())
        .def("insert",static_cast<void (Array1D<double>::*)(Array1D<double>&,int)>(&Array1D<double>::insert))
        .def("insert",static_cast<void (Array1D<double>::*)(const double&,int)>(&Array1D<double>::insert))
        .def("erase",&Array1D<double>::erase)
        .def("DumpBinary",static_cast<void (Array1D<double>::*)(FILE*) const>(&Array1D<double>::DumpBinary))
        .def("DumpBinary",static_cast<void (Array1D<double>::*)(char*)>(&Array1D<double>::DumpBinary))
        .def("ReadBinary",static_cast<void (Array1D<double>::*)(FILE*)>(&Array1D<double>::ReadBinary))
        .def("ReadBinary",static_cast<void (Array1D<double>::*)(char*)>(&Array1D<double>::ReadBinary))
        .def("pyElement",&Array1D<double>::operator[])
        .def("ReadBinary4py",&Array1D<double>::ReadBinary4py)
        .def("DumpBinary4py",&Array1D<double>::DumpBinary4py)
        .def("setArray",&Array1D<double>::setArray)
        .def("setnpdblArray",py::vectorize(&Array1D<double>::setnpdblArray))
        //.def("getnpdblArray",py::vectorize(&Array1D<double>::getnpdblArray))
        .def("flatten",&Array1D<double>::flatten)
        .def("type",&Array1D<double>::type)
        .def("shape",&Array1D<double>::shape)
        .def("assign",&Array1D<double>::assign)
        ;

      py::class_<Array2D<int>>(m,"intArray2D")
        .def(py::init<>())
        .def(py::init<const int&,const int&>())
        .def(py::init<const int&,const int&,const int&>())
        .def(py::init<const Array2D<int> &>())
        .def("Clear",&Array2D<int>::Clear)
        .def("XSize",&Array2D<int>::XSize)
        .def("YSize",&Array2D<int>::YSize)
        .def("Resize",static_cast<void (Array2D<int>::*)(const int&,const int&)>(&Array2D<int>::Resize))
        .def("Resize",static_cast<void (Array2D<int>::*)(const int&,const int&,const int&)>(&Array2D<int>::Resize))
        .def("SetValue",&Array2D<int>::SetValue)
        .def("GetArrayPointer",&Array2D<int>::GetArrayPointer)
        .def("GetConstArrayPointer",&Array2D<int>::GetConstArrayPointer)
        .def("__getitem__", [](Array2D<int> a, int b,int c) {
          return a(b,c);
        }, py::is_operator())
        .def("insertRow",static_cast<void (Array2D<int>::*)(Array1D<int>&,int)>(&Array2D<int>::insertRow))
        .def("insertRow",static_cast<void (Array2D<int>::*)(Array2D<int>&,int)>(&Array2D<int>::insertRow))
        .def("eraseRow",&Array2D<int>::eraseRow)
        .def("insertCol",static_cast<void (Array2D<int>::*)(Array1D<int>&,int)>(&Array2D<int>::insertCol))
        .def("insertCol",static_cast<void (Array2D<int>::*)(Array2D<int>&,int)>(&Array2D<int>::insertCol))
        .def("eraseCol",&Array2D<int>::eraseCol)
        .def("DumpBinary",static_cast<void (Array2D<int>::*)(FILE*) const>(&Array2D<int>::DumpBinary))
        .def("DumpBinary",static_cast<void (Array2D<int>::*)(char*)>(&Array2D<int>::DumpBinary))
        .def("ReadBinary",static_cast<void (Array2D<int>::*)(FILE*)>(&Array2D<int>::ReadBinary))
        .def("ReadBinary",static_cast<void (Array2D<int>::*)(char*)>(&Array2D<int>::ReadBinary))
        .def("pyElement",&Array2D<int>::operator[])
        .def("getRow",&Array2D<int>::getRow)
        .def("ReadBinary4py",&Array2D<int>::ReadBinary4py)
        .def("DumpBinary4py",&Array2D<int>::DumpBinary4py)
        .def("setArray",&Array2D<int>::setArray)
        .def("flatten",&Array2D<int>::flatten)
        .def("type",&Array2D<int>::type)
        //.def("setnpintArray",py::vectorize(&Array2D<int>::setnpintArray))
        //.def("getnpintArray",py::vectorize(&Array2D<int>::getnpintArray))
        .def("shape",&Array2D<int>::shape)
        .def("at",&Array2D<int>::at)
        .def("assign",&Array2D<int>::assign)
        ;

      py::class_<Array2D<double>>(m,"dblArray2D")
        .def(py::init<>())
        .def(py::init<const int&,const int&>())
        .def(py::init<const int&,const int&,const int&>())
        .def(py::init<const Array2D<double> &>())
        .def("Clear",&Array2D<double>::Clear)
        .def("XSize",&Array2D<double>::XSize)
        .def("YSize",&Array2D<double>::YSize)
        .def("Resize",static_cast<void (Array2D<double>::*)(const int&,const int&)>(&Array2D<double>::Resize))
        .def("Resize",static_cast<void (Array2D<double>::*)(const int&,const int&,const double&)>(&Array2D<double>::Resize))
        .def("SetValue",&Array2D<double>::SetValue)
        .def("GetArrayPointer",&Array2D<double>::GetArrayPointer)
        .def("GetConstArrayPointer",&Array2D<double>::GetConstArrayPointer)
        .def("__getitem__", [](Array2D<double> a, int b,int c) {
          return a(b,c);
        }, py::is_operator())
        .def("insertRow",static_cast<void (Array2D<double>::*)(Array1D<double>&,int)>(&Array2D<double>::insertRow))
        .def("insertRow",static_cast<void (Array2D<double>::*)(Array2D<double>&,int)>(&Array2D<double>::insertRow))
        .def("eraseRow",&Array2D<double>::eraseRow)
        .def("insertCol",static_cast<void (Array2D<double>::*)(Array1D<double>&,int)>(&Array2D<double>::insertCol))
        .def("insertCol",static_cast<void (Array2D<double>::*)(Array2D<double>&,int)>(&Array2D<double>::insertCol))
        .def("eraseCol",&Array2D<double>::eraseCol)
        .def("DumpBinary",static_cast<void (Array2D<double>::*)(FILE*) const>(&Array2D<double>::DumpBinary))
        .def("DumpBinary",static_cast<void (Array2D<double>::*)(char*)>(&Array2D<double>::DumpBinary))
        .def("ReadBinary",static_cast<void (Array2D<double>::*)(FILE*)>(&Array2D<double>::ReadBinary))
        .def("ReadBinary",static_cast<void (Array2D<double>::*)(char*)>(&Array2D<double>::ReadBinary))
        .def("pyElement",&Array2D<double>::operator[])
        .def("getRow",&Array2D<double>::getRow)
        .def("ReadBinary4py",&Array2D<double>::ReadBinary4py)
        .def("DumpBinary4py",&Array2D<double>::DumpBinary4py)
        .def("setArray",&Array2D<double>::setArray)
        .def("flatten",&Array2D<double>::flatten)
        .def("type",&Array2D<double>::type)
        //.def("setnpdblArray",py::vectorize(&Array2D<double>::setnpdblArray))
        //.def("getnpdblArray",py::vectorize(&Array2D<double>::getnpdblArray))
        .def("shape",&Array2D<double>::shape)
        .def("at",&Array2D<double>::at)
        .def("assign",&Array2D<double>::assign)
        ;

      m.def("read_datafile",&read_datafile<int>);
      m.def("read_datafile",&read_datafile<double>);
      m.def("read_datafile_1d",&read_datafile_1d<int>);
      m.def("read_datafile_1d",&read_datafile_1d<double>);
      m.def("write_datafile_size",&write_datafile_size<int>);
      m.def("write_datafile_size",&write_datafile_size<double>);
      m.def("write_datafile_1d",&write_datafile_1d<int>);
      m.def("write_datafile_1d",&write_datafile_1d<double>);
      m.def("read_datafileVS",&foo<int>);
      m.def("read_datafileVS",&foo<double>);
      m.def("read_datafileVS",&fun<double>);
      m.def("read_datafileVS",&fun<int>);
      m.def("read_datafileVS",&dummy<int>);
      m.def("read_datafileVS",&dummy<double>);
      m.def("write_datafile",&rush<int>);
      m.def("write_datafile",&rush<double>);
      m.def("write_datafile",&stylebender<int>);
      m.def("write_datafile",&stylebender<double>);
      m.def("write_datafile",&rowdy<int>);
      m.def("write_datafile",&rowdy<double>);

      m.def("array1Dto2D",&array1Dto2D<int>);
      m.def("array1Dto2D",&array1Dto2D<double>);
      m.def("array2Dto1D",&array2Dto1D<int>);
      m.def("array2Dto1D",&array2Dto1D<double>);
      m.def("generate_multigrid",&generate_multigrid<int>);
      m.def("generate_multigrid",&generate_multigrid<double>);
      m.def("transpose",&transpose<double>);
      m.def("transpose",&transpose<int>);
      m.def("flatten",&flatten);
      m.def("fold_1dto2d_rowfirst",&fold_1dto2d_rowfirst);
      m.def("fold_1dto2d_colfirst",&fold_1dto2d_colfirst);
      m.def("access",&accessPythonHelper);
      m.def("getRow",&getRow<double>);
      m.def("getRow",&getRow<int>);
      m.def("getCol",&getCol<int>);
      m.def("getCol",&getCol<double>);
      m.def("subVector",&subVector<int>);
      m.def("subVector",&subVector<double>);
      m.def("subMatrix_row_int",&subMatrix_row<int>);
      m.def("subMatrix_row_dbl",&subMatrix_row<double>);
      m.def("subMatrix_col",&subMatrix_col<double>);
      m.def("subMatrix_col",&subMatrix_col<int>);
      m.def("matPvec",&matPvec<double>);
      m.def("matPvec",&matPvec<int>);
      m.def("maxVal",&maxVal<double>);
      m.def("maxVal",&maxVal<int>);
      m.def("setdiff",&setdiff);
      m.def("setdiff_s",&setdiff_s);
      m.def("shell_sort_col",&shell_sort_col);
      m.def("shell_sort_all",&shell_sort_all);
      m.def("find",&find<double>);
      m.def("find",&find<int>);
      m.def("prodAlphaMatVec",&prodAlphaMatVec);
      m.def("prodAlphaMatTVec",&prodAlphaMatTVec);
      m.def("prodAlphaMatMat",&prodAlphaMatMat);
      m.def("prodAlphaMatTMat",&prodAlphaMatTMat);
      m.def("addVecAlphaVecPow",&addVecAlphaVecPow);
      m.def("prod_vecTmatvec",&prod_vecTmatvec);
      m.def("MatTMat",&MatTMat);
      m.def("delRow",&delRow<double>);
      m.def("delRow",&delRow<int>);
      m.def("paddMatColScal",&paddMatColScal);
      m.def("vecIsInArray",&vecIsInArray);
      m.def("select_kth",&select_kth);
      m.def("logdeterm",&logdeterm);
      m.def("trace",&trace);
      m.def("evalLogMVN",&evalLogMVN);
      m.def("diag",&diag);
      m.def("mtxdel",&mtxdel);
      m.def("norm",static_cast<double (*)(Array1D<double> &)>(&norm));
      m.def("dotT",&dotT);
      m.def("INV",&INV);
      m.def("AinvH",&AinvH);
      m.def("Ainvb",&Ainvb);
      m.def("LSTSQ",&LSTSQ);
      m.def("QR",&QR);
      m.def("SVD",&SVD);
      m.def("dist_sq",&dist_sq);
      m.def("Trans",&Trans);
      m.def("paste",static_cast<void (*)(Array1D<double> &, Array1D<double> &, Array2D<double> &)>(&paste));
      m.def("paste",static_cast<void (*)(Array1D<int> &, Array1D<int> &, Array2D<int> &)>(&paste));
      m.def("paste",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&)>(&paste));
      m.def("merge",static_cast<void (*)(Array1D<double>&, Array1D<double>&, Array1D<double>&)>(&merge));
      m.def("merge",static_cast<void (*)(Array2D<double>&, Array2D<double>&, Array2D<double>&)>(&merge));
      m.def("merge",static_cast<void (*)(Array1D<int>&, Array1D<int>&, Array1D<int>&)>(&merge));
      m.def("append",static_cast<void (*)(Array1D<double>&, Array1D<double>&)>(&append));
      m.def("append",static_cast<void (*)(Array1D<int>&, Array1D<int>&)>(&append));
      m.def("swap",static_cast<void (*)(Array1D<double>&,int,int)>(&swap));
      m.def("swap",static_cast<void (*)(Array2D<double>&,int,int)>(&swap));
      m.def("addVal",static_cast<void (*)(int, double *, double)>(&addVal));
      m.def("addVal",static_cast<void (*)(int, int *, int)>(&addVal));
      m.def("addVal",static_cast<void (*)(Array1D<double> &, double )>(&addVal));
      m.def("addVal",static_cast<void (*)(Array1D<int> &, int)>(&addVal));
      m.def("addVal",static_cast<void (*)(Array2D<double> &, double )>(&addVal));
      m.def("addVal",static_cast<void (*)(Array2D<int> &, int)>(&addVal));
      m.def("shell_sort",static_cast<void (*)(int *, int)>(&shell_sort));
      m.def("shell_sort",static_cast<void (*)(Array1D<int>&)>(&shell_sort));
      m.def("shell_sort",static_cast<void (*)(Array1D<double>&)>(&shell_sort));
      m.def("quicksort3",static_cast<void (*)(Array1D<double>&,int,int)>(&quicksort3));
      m.def("quicksort3",static_cast<void (*)(Array2D<double>&,int,int,int)>(&quicksort3));
      m.def("quicksort3",static_cast<void (*)(Array2D<double>&,int,int)>(&quicksort3));
      m.def("intersect",static_cast<void (*)(Array1D<int> &, Array1D<int> &, Array1D<int> &,Array1D<int> &,Array1D<int> &)>(&intersect));
      m.def("intersect",static_cast<void (*)(Array1D<int> &, Array1D<int> &, Array1D<int> &)>(&intersect));
      m.def("delCol",static_cast<void (*)(Array2D<double>&,int)>(&delCol));
      m.def("delCol",static_cast<void (*)(Array2D<int>&,int)>(&delCol));
      m.def("delCol",static_cast<void (*)(Array1D<double>&,int)>(&delCol));
      m.def("delCol",static_cast<void (*)(Array1D<int>&,int)>(&delCol));
      m.def("paddMatRow",static_cast<void (*)(Array2D<double>&, Array1D<double>&)>(&paddMatRow));
      m.def("paddMatRow",static_cast<void (*)(Array2D<int>&, Array1D<int>&)>(&paddMatRow));
      m.def("paddMatCol",static_cast<void (*)(Array2D<double>&, Array1D<double>&)>(&paddMatCol));
      m.def("paddMatCol",static_cast<void (*)(Array2D<int>&, Array1D<int>&)>(&paddMatCol));
      m.def("is_equal",static_cast<bool (*)(Array1D<int>&,Array1D<int>&)>(&is_equal));
      m.def("is_equal",static_cast<bool (*)(Array1D<double>&,Array1D<double>&)>(&is_equal));
      m.def("is_less",static_cast<bool (*)(Array1D<double>&,Array1D<double>&)>(&is_less));
      m.def("is_less",static_cast<bool (*)(Array1D<int>&,Array1D<int>&)>(&is_less));
      m.def("copy",static_cast<Array1D<double> (*)(Array1D<double>&)>(&copy));
      m.def("copy",static_cast<Array2D<double> (*)(Array2D<double>&)>(&copy));
      m.def("add",static_cast<Array1D<double> (*)(Array1D<double>&,Array1D<double>&)>(&add));
      m.def("add",static_cast<Array2D<double> (*)(Array2D<double>&,Array2D<double>&)>(&add));
      m.def("addinplace",static_cast<void (*)(Array2D<double>&,Array2D<double>&)>(&addinplace));
      m.def("addinplace",static_cast<void (*)(Array1D<double>&,Array1D<double>&)>(&addinplace));
      m.def("subtract",static_cast<Array1D<double> (*)(Array1D<double>&,Array1D<double>&)>(&subtract));
      m.def("subtract",static_cast<Array2D<double> (*)(Array2D<double>&,Array2D<double>&)>(&subtract));
      m.def("subtractinplace",static_cast<void (*)(Array2D<double>&,Array2D<double>&)>(&subtractinplace));
      m.def("subtractinplace",static_cast<void (*)(Array1D<double>&,Array1D<double>&)>(&subtractinplace));
      m.def("scale",static_cast<Array1D<double> (*)(Array1D<double>&,double)>(&scale));
      m.def("scale",static_cast<Array2D<double> (*)(Array2D<double>&,double)>(&scale));
      m.def("scaleinplace",static_cast<void (*)(Array1D<double>&,double)>(&scaleinplace));
      m.def("scaleinplace",static_cast<void (*)(Array1D<int>&,int)>(&scaleinplace));
      m.def("scaleinplace",static_cast<void (*)(Array2D<double>&, double)>(&scaleinplace));
      m.def("scaleinplace",static_cast<void (*)(Array2D<int>&,int)>(&scaleinplace));
      m.def("dotdivide",static_cast<Array2D<double> (*)(Array2D<double>&,Array2D<double>&)>(&dotdivide));
      m.def("dotdivide",static_cast<Array1D<double> (*)(Array1D<double>&,Array1D<double>&)>(&dotdivide));
      m.def("dotmult",static_cast<Array2D<double> (*)(Array2D<double>&,Array2D<double>&)>(&dotmult));
      m.def("dotmult",static_cast<Array1D<double> (*)(Array1D<double>&,Array1D<double>&)>(&dotmult));
      m.def("dot",static_cast<double (*)(Array1D<double>&,Array1D<double>&)>(&dot));
      m.def("dot",&method);
      m.def("dot",&swag);
      m.def("dot",static_cast<Array2D<double> (*)(Array2D<double>&, Array2D<double>&)>(&dot));
      m.def("printarray",static_cast<void (*)(Array1D<double>&)>(&printarray));
      m.def("printarray",static_cast<void (*)(Array1D<int>&)>(&printarray));
      m.def("printarray",static_cast<void (*)(Array2D<double>&)>(&printarray));
      m.def("printarray",static_cast<void (*)(Array2D<int>&)>(&printarray));

      m.def("setnpdblArray",&setnpdblArray);
      m.def("setnpintArray",&setnpintArray);
      m.def("getnpdblArray",static_cast<py::array_t<double> (*)(Array1D<double> &)>(&getnpdblArray),py::arg("x"),py::return_value_policy::take_ownership);
      m.def("getnpintArray",static_cast<py::array_t<int> (*)(Array1D<int> &)>(&getnpintArray),py::arg("x"),py::return_value_policy::take_ownership);
      m.def("getnpdblArray",static_cast<vector<double> (*)(Array2D<double> &)>(&getnpdblArray),py::arg("x"),py::return_value_policy::take_ownership);
      m.def("getnpintArray",static_cast<vector<int> (*)(Array2D<int> &)>(&getnpintArray),py::arg("x"),py::return_value_policy::take_ownership);




}
