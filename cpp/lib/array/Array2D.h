/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file Array2D.h
/// \brief 2D Array class for any type T


#ifndef ARRAY2D_H_SEEN
#define ARRAY2D_H_SEEN

#include <stddef.h>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <typeinfo>
#include "Array1D.h"

using namespace std;

/// \class  Array2D
/// \brief  Stores data of any type T in a 2D array
///
/// This class also provides a Fortran-like access operator ()
/// as well as a function to access the data in the array through a pointer that
/// can be passed to F77 or C routines.
/// \author  Bert Debusschere <bjdebus@sandia.gov>
/// \date  Jan 2005
/// \note  Inspired by Helgi Adalsteinsson's Array class implementation
/// \todo  Define copy constructor
// COLUMN MAJOR FORMAT

template<typename T>
class Array2D{
private:

public:
  // These two quantities used to be private but making them public
  // allows for easy access to python interface as a "list"
  int xsize_;
  int ysize_;
  vector<T> data_;
  Array1D<T> arraycopy;
  Array1D<T> rowvec;

  /// \brief Default constructor, which does not allocate any memory
  Array2D(): xsize_(0), ysize_(0) {};

  /// \brief Constructor that allocates the memory
  Array2D(const int& nx, const int& ny):  xsize_(nx), ysize_(ny){
  data_.resize(xsize_*ysize_);
  }

  /// \brief Constructor that allocates and initializes the data to a constant t
  Array2D(const int& nx, const int& ny, const T& t):  xsize_(nx), ysize_(ny){
    data_.resize(xsize_*ysize_ , t);
  }

  /// \brief Copy constructor
  Array2D(const Array2D &obj): xsize_(obj.xsize_), ysize_(obj.ysize_), data_(obj.data_) {};

  /// \brief Destructor that frees up the memory
  ~Array2D() {data_.clear();}

  /// \brief Function to clear the memory
  void Clear() {
  xsize_ = 0;
  ysize_ = 0;
  data_.clear();
  }

  /// \brief Returns size in the x-direction
  int XSize() const {return xsize_;}
  /// \brief Returns size in the y-direction
  int YSize() const {return ysize_;}

  /// \brief Resizes the array
  /// \warning In its current implementation, most of the original data
  void Resize(const int& nx, const int& ny) {
    xsize_ = nx;
    ysize_ = ny;
    data_.resize(xsize_*ysize_);
  }

  /// \brief Resizes the array and sets ALL entries to the specified value
  /// \warning All original data will get lost if this function is used!
  void Resize(const int& nx, const int& ny, const T& t) {
  data_.clear();
    xsize_ = nx;
    ysize_ = ny;
    data_.resize(xsize_*ysize_, t);
  }

  /// \brief Set all values in the array to the given value
  void SetValue(const T& t){
  for(int i=0; i < data_.size(); i++){
    data_[i] = t;
  }
  }

  /// \brief Return a pointer to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  T* GetArrayPointer() {
    return &(data_[0]);
  }

  /// \brief Return a cont point to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  const T* GetConstArrayPointer() const {
    return &(data_[0]);
  }

  /// \brief C-like () operator to access values in the 2D data array
  // values accessed in a row-major format
  T& operator()(int ix,int iy) {return data_[ix + xsize_*iy];}
  const T& operator()(int ix,int iy) const {return data_[ix + xsize_*iy];}

   /// \brief Insert array insarr as a row into position ix
  void insertRow(Array1D<T>& insarr,int ix){
    if (ix<0 || ix>xsize_)
      throw Tantrum("Array2D:insertRow():: insert index out of bounds.");
    if ( insarr.Length() != ysize_ )
      throw Tantrum("Array2D:insertRow():: insert row size does not match.");

    vector<T> data_old;
    data_old=data_;

    xsize_ += 1; // new number of rows
    data_.resize(xsize_*ysize_);

    for(int iy=0;iy<ysize_;iy++){
      for(int i=0; i < ix; i++)
        data_[i+xsize_*iy] = data_old[i+(xsize_-1)*iy];
      data_[ix+xsize_*iy]=insarr(iy);
      for(int i=ix+1; i < xsize_; i++)
        data_[i+xsize_*iy] = data_old[i-1+(xsize_-1)*iy];
    }
  }

  /// \brief Insert a 2d-array insarr into a row position ix
  void insertRow(Array2D<T>& insarr,int ix){
    if (ix<0 || ix>xsize_)
      throw Tantrum("Array2D:insertRow():: insert index out of bounds.");
    if ( insarr.YSize() != ysize_ )
      throw Tantrum("Array2D:insertRow():: insert row size does not match.");

    vector<T> data_old;
    data_old=data_;

    int insx=insarr.XSize();

    xsize_ += insx;
    data_.resize(xsize_*ysize_);

    for(int iy=0;iy<ysize_;iy++){
      for(int i=0; i < ix; i++)
      data_[i+xsize_*iy] = data_old[i+(xsize_-insx)*iy];
      for(int i=ix; i < ix+insx; i++)
      data_[i+xsize_*iy]=insarr(i-ix,iy);
      for(int i=ix+insx; i < xsize_; i++)
      data_[i+xsize_*iy] = data_old[i-insx+(xsize_-insx)*iy];
    }
  }

  /// \brief Erase the row ix
  void eraseRow(int ix){
    if (ix<0 || ix>=xsize_)
      throw Tantrum("Array2D:eraseRow():: erase index out of bounds.");

    vector<T> data_old;
    data_old=data_;

    xsize_-=1;
    data_.resize(xsize_*ysize_);

    for(int iy=0;iy<ysize_;iy++){
      for(int i=0; i < ix; i++)
      data_[i+xsize_*iy] = data_old[i+(xsize_+1)*iy];
      for(int i=ix; i < xsize_; i++)
      data_[i+xsize_*iy] = data_old[i+1+(xsize_+1)*iy];
    }

    //if (xsize_==0)
    //  printf("eraseRow(): WARNING: the xsize is zeroed!");

  }

  // /// \brief Insert array insarr as a column into position iy
  void insertCol(Array1D<T>& insarr,int iy){
    if (iy<0 || iy>ysize_)
      throw Tantrum("Array2D:insertCol():: insert index out of bounds.");
    if ( insarr.Length() != xsize_ )
      throw Tantrum("Array2D:insertCol():: insert column size does not match.");


    T* ptr=insarr.GetArrayPointer();
    data_.insert(data_.begin()+xsize_*iy,ptr,ptr+xsize_);

    ysize_+=1;

  }

  /// \brief Insert a 2d-array insarr into a column position iy
  void insertCol(Array2D<T>& insarr,int iy){
    if (iy<0 || iy>ysize_)
      throw Tantrum("Array2D:insertCol():: insert index out of bounds.");
    if ( insarr.XSize() != xsize_ )
      throw Tantrum("Array2D:insertRow():: insert column size does not match.");

    int insy=insarr.YSize();

    T* ptr=insarr.GetArrayPointer();
    data_.insert(data_.begin()+xsize_*iy,ptr,ptr+xsize_*insy);

    ysize_+=insy;
    }

    /// \brief Erase the column iy
  void eraseCol(int iy){
    if (iy<0 || iy>=ysize_)
      throw Tantrum("Array2D:eraseCol():: erase index out of bounds.");

    data_.erase(data_.begin()+xsize_*iy,data_.begin()+xsize_*(iy+1));

    ysize_-=1;

    //if (ysize_==0)
    // printf("eraseCol(): WARNING: the ysize is zeroed!");

  }

   /// \brief Dump contents of the array to a file in binary format
  void DumpBinary(FILE* f_out) const {
    fwrite(&xsize_,sizeof(xsize_),1,f_out);
    fwrite(&ysize_,sizeof(ysize_),1,f_out);
    fwrite(this->GetConstArrayPointer(),sizeof(T),xsize_*ysize_,f_out);
  }


  /// \brief Read contents of the array from a file in binary format
  void ReadBinary(FILE* f_in){
    fread(&xsize_,sizeof(xsize_),1,f_in);
    fread(&ysize_,sizeof(ysize_),1,f_in);
    data_.resize(xsize_*ysize_);
    fread(this->GetArrayPointer(),sizeof(T),xsize_*ysize_,f_in);
  }

  /********************************************************
  // Methods for interfacing with python
  ********************************************************/

  // assignment operator []
  // allows for calling Array2D using [i][j] notation
  // make more efficient by setting two vectors equal
  Array1D<T>& operator[](int ix) {
    // get the ith row
    int stride = xsize_;
    rowvec.Resize(ysize_);
    for (int iy = 0; iy < ysize_; iy++){
      rowvec(iy) = data_[ix + stride*iy];
    }
    return rowvec;
  }

  // For calling shape in Python
  vector<int> shape(){
    vector<int> s (2,0);
    s[0] = this -> XSize();
    s[1] = this -> YSize();
    return s;
  }

  void getRow(int row){
    arraycopy.Resize(ysize_,0);
    int stride = xsize_;
    for (int i = 0; i < ysize_; i++){
      arraycopy[i] = data_[i*stride + row];
    }
  }

  // read binary file created with DumpBinary
  // Cannot use numpy's from files
  // only for use in c++
  void DumpBinary(char *filename){
  FILE *f_out;
  f_out = fopen(filename,"wb");
  fwrite(&xsize_,sizeof(xsize_),1,f_out);
  fwrite(&ysize_,sizeof(ysize_),1,f_out);
  fwrite(this->GetConstArrayPointer(),sizeof(T),xsize_*ysize_,f_out);
  fclose(f_out);
  }

  // Only for use if DumpBinary was used
  // can only be read in c++
  // can be opened with ReadBinary(FILE* file) above
  void ReadBinary(char *filename){
  FILE *f_in;
  f_in = fopen(filename,"rb");
  fread(&xsize_,sizeof(xsize_),1,f_in);
  fread(&ysize_,sizeof(ysize_),1,f_in);
  data_.resize(xsize_*ysize_);
  fread(this->GetArrayPointer(),sizeof(T),xsize_*ysize_,f_in);
  fclose(f_in);
  }

  // creates binary file that can be read with numpy's fromfile
  void DumpBinary4py(char *filename){
  ofstream f_out;
  f_out.open(filename, ios::out | ios::binary);
  f_out.write((char*)this->GetArrayPointer(),sizeof(T[xsize_*ysize_])); // convert array pointer to char string
  f_out.close();
  }

  // can read in DumpBinary4py output, but needs size of vector
  // fromfile can automatically detect size file, so, if need by, one can use numpy's fromfile to determine # of elements
  void ReadBinary4py(char *filename, int n1, int n2){
  xsize_ = n1;
  ysize_ = n2;
  ifstream f_in;
  f_in.open(filename, ios::in | ios::binary);
  f_in.read((char*)this->GetArrayPointer(),sizeof(T[xsize_*ysize_])); // convert array pointer to char string
  f_in.close();
  }

  // Set user-defined list to data_ vector
  // This will work even for string type
  void setArray(vector<T> inarray){
    data_ = inarray;
    // xsize_ = inarray.size();
  }

  // Sets user-defined 2d numpy array to data_ vector
  // This is not to be used for a string type
  void setnpdblArray(double* inarray, int n1, int n2){
    xsize_ = n1;
    ysize_ = n2;
    data_.assign(inarray,inarray+n1*n2);
  }

  // get numpy double array from data_ vector
  void getnpdblArray(double* outarray){
    copy(data_.begin(), data_.end(), outarray);
  }

  // Sets user-defined 2d numpy array to data_ vector
  // This is not to be used for a string type
  void setnpintArray(long* inarray, int n1, int n2){
    xsize_ = n1;
    ysize_ = n2;
    data_.assign(inarray,inarray+n1*n2);
  }

  // get numpy double array from data_ vector
  void getnpintArray(long* outarray){
    copy(data_.begin(), data_.end(), outarray);
  }

  // Get the value at location x,y
  T& at(int ix,int iy){
    return data_[ix + xsize_*iy];
  }

  // Returns data_ vector as a list in python in row-major (?)
  // Also acts as a print to see individual elements
  vector<T> flatten(){
    return data_;
  }

  string type(){
    const char* s = typeid(data_[0]).name();
    if (string(s) == string("Ss") ){
      return "string";
    }
    else if (strcmp(s,"i") == 0){
      return "int";
    }
    else {
      return "double";
    }
  }

  // For python, allows the user to assign a value to a specific index (x,y)
  void assign(const int x,const int y,const T val){
    data_[x + xsize_*y] = val;
  }
};

#endif /* ARRAY2D_H_SEEN */
