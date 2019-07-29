/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file Array1D.h
/// \brief 1D Array class for any type T

#ifndef ARRAY1D_H_SEEN
#define ARRAY1D_H_SEEN

#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <typeinfo>

#include "error_handlers.h"

using namespace std;

// template<typename T> T max_test(T a, T b) { return a > b ? a : b; }

/// \class  Array1D
/// \brief  Stores data of any type T in a 1D array
///
/// This class also provides a Fortran-like access operator ()
/// as well as a function to access the data in the array through a pointer that
/// can be passed to F77 or C routines.
/// \author Bert Debusschere <bjdebus@sandia.gov>
/// \date Apr 2005 - Nov 2007
/// \note Inspired by Helgi Adalsteinsson's Array class implementation
/// \todo double check copy constructor
//column major for fortran blas
template<typename T>
class Array1D{
private:

public:
  // These two quantities used to be private but making them public
  // allows for easy access to python interface as a "list"
  int xsize_; // public (used to be private) size of vector
  vector<T> data_; // public (used to be private) copy of data vector

  /// Default constructor, which does not allocate any memory
  Array1D(): xsize_(0) {};

  /// \brief Constructor that allocates the memory
  Array1D(const int& nx): xsize_(nx) {
    data_.resize(xsize_);
  }

  /// Constructor that allocates and initializes the data to a value t
  Array1D(const int& nx, const T& t): xsize_(nx) {
    data_.resize(xsize_, t);
  }

  /// \brief Assignment operator copies the data structure by value
  Array1D& operator=(const Array1D &obj) {
    xsize_ = obj.xsize_;
    data_ = obj.data_;
    return *this;
  }

  /// \brief Copy constructor
  Array1D(const Array1D &obj): xsize_(obj.xsize_), data_(obj.data_) {};

  /// Destructor that frees up the memory
  ~Array1D() {data_.clear();}

  /// \brief Function to clear the memory
  void Clear() {
    xsize_ = 0;
    data_.clear();
  }

  /// \brief Returns size in the x-direction
  int XSize() const {return xsize_;}

  /// Returns length (i.e. size in the x-direction)
  int Length() const {return xsize_;}

  /// \brief Resizes the array
  void Resize(const int& nx) {
    xsize_ = nx;
    data_.resize(xsize_);
  }

  /// \brief Resizes the array and sets ALL entries to the specified value
  /// \warning All original data will get lost if this function is used!
  void Resize(const int& nx, const T& t) {
    data_.clear();
    xsize_ = nx;
    data_.resize(xsize_, t);
  }

  /// \brief Set all values in the array to the given value
  void SetValue(const T& t){
    for(int i=0; i < data_.size(); i++){
      data_[i] = t;
    }
  }

  /// \brief Add element to the end of the vector
  void PushBack(const T& t){
     xsize_ += 1;
     data_.push_back(t);
  }

  /// \brief Return a pointer to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  T* GetArrayPointer() {
    return &(data_[0]);
  }

  /// \brief Return a const point to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  const T* GetConstArrayPointer() const {
    return &(data_[0]);
  }

  // allows access element by element, e.g. this(i) gives data_[i]
  T& operator()(int ix) {return data_[ix];}
  const T& operator()(int ix) const {return data_[ix];}

  /// \brief Insert a given array to the position ix
  /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
  void insert(Array1D<T>& insarr,int ix){
    if (ix<0 || ix>xsize_){
      throw Tantrum("Array1D:insert():: insert index out of bounds.");
    }
    int addsize = insarr.Length();
    xsize_+=addsize;
    T* ptr=insarr.GetArrayPointer();
    data_.insert(data_.begin()+ix,ptr,ptr+addsize);
  }

  /// \brief Insert a given value to the position ix
  /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
  void insert(const T& insval,int ix){
    if (ix<0 || ix>xsize_)
      throw Tantrum("Array1D:insert():: insert index out of bounds.");
    xsize_+=1;
    data_.insert(data_.begin()+ix,insval);
  }

  /// \brief Erase the value from the position ix
  void erase(int ix){
    if (ix<0 || ix>=xsize_)
      throw Tantrum("Array1D:erase():: erase index out of bounds.");
    xsize_-=1;
    data_.erase(data_.begin()+ix);
  }

  /// \brief Dump contents of the array to a file in binary format
  void DumpBinary(FILE* f_out) const {
    fwrite(&xsize_,sizeof(xsize_),1,f_out);
    fwrite(this->GetConstArrayPointer(),sizeof(T),xsize_,f_out);
  }

  /// \brief Read contents of the array from a file in binary format
  void ReadBinary(FILE* f_in){
    fread(&xsize_,sizeof(xsize_),1,f_in);
    data_.resize(xsize_);
    fread(this->GetArrayPointer(),sizeof(T),xsize_,f_in);
  }

  /**********************************************************
  // Methods for interfacing with python 
  **********************************************************/

  // For calling [] in Python
  T& operator[](int i) {return data_[i];}

  /// \brief Dump contents of the array to a file in binary format
  // cannot be read with numpy's fromfile after creation
  void DumpBinary(char *filename){
    FILE *f_out; 
    f_out = fopen(filename,"wb");
    fwrite(&xsize_,sizeof(xsize_),1,f_out);
    fwrite(this->GetConstArrayPointer(),sizeof(T),xsize_,f_out);
    fclose(f_out);
  }

  // read binary file created with DumpBinary
  // Cannot use numpy's from files
  // only for use in c++
  void ReadBinary(char *filename){
    FILE *f_in;
    f_in = fopen(filename,"rb");
    fread(&xsize_,sizeof(xsize_),1,f_in);
    data_.resize(xsize_);
    fread(this->GetArrayPointer(),sizeof(T),xsize_,f_in);
    fclose(f_in);
  }

  // Following two methods are not compatable with certain clang comilers
  // creates binary file that can be read with numpy's fromfile 
  void DumpBinary4py(char *filename){
    ofstream f_out;
    f_out.open(filename, ios::out | ios::binary);
    f_out.write((char*)this->GetArrayPointer(),sizeof(T[xsize_])); // convert array pointer to char string
    f_out.close();
  }

  // can read in DumpBinary4py output, but needs size of vector
  // fromfile can automatically detect size file, so, if need by, one can use numpy's fromfile to determine # of elements
  void ReadBinary4py(char *filename, int n){
    xsize_ = n; 
    ifstream f_in;
    f_in.open(filename, ios::in | ios::binary);
    f_in.read((char*)this->GetArrayPointer(),sizeof(T[xsize_])); // convert array pointer to char string
    f_in.close();
  }

  // Set user-defined list to data_ vector
  // This will work even for string type
  void setArray(vector<T> inarray){
    data_ = inarray; 
    xsize_ = inarray.size();
  }

  // Returns data_ vector as a list in python
  // Also acts as a print to see individual elements
  vector<T> flatten(){
    return data_;
  }

  string type(){
    return "string";
  }
};

template<>
class Array1D <int>{
private:
  int xsize_; // private size of vector
public:
  vector<int> data_; // private copy of data vector

  /// Default constructor, which does not allocate any memory
  Array1D(): xsize_(0) {};

  /// \brief Constructor that allocates the memory
  Array1D(const int& nx): xsize_(nx) {
    data_.resize(xsize_);
  }

  /// Constructor that allocates and initializes the data to a value t
  Array1D(const int& nx, const int& t): xsize_(nx) {
    data_.resize(xsize_, t);
  }

  /// \brief Assignment operator copies the data structure by value
  Array1D& operator=(const Array1D &obj) {
    xsize_ = obj.xsize_;
    data_ = obj.data_;
    return *this;
  }

  /// \brief Copy constructor
  Array1D(const Array1D &obj): xsize_(obj.xsize_), data_(obj.data_) {};

  /// Destructor that frees up the memory
  ~Array1D() {data_.clear();}

  /// \brief Function to clear the memory
  void Clear() {
    xsize_ = 0;
    data_.clear();
  }

  /// \brief Returns size in the x-direction
  int XSize() const {return xsize_;}

  /// Returns length (i.e. size in the x-direction)
  int Length() const {return xsize_;}

  /// \brief Resizes the array
  void Resize(const int& nx) {
    xsize_ = nx;
    data_.resize(xsize_);
  }

  /// \brief Resizes the array and sets ALL entries to the specified value
  /// \warning All original data will get lost if this function is used!
  void Resize(const int& nx, const int& t) {
    data_.clear();
    xsize_ = nx;
    data_.resize(xsize_, t);
  }

  /// \brief Set all values in the array to the given value
  void SetValue(const int& t){
    for(int i=0; i < (int)data_.size(); i++){
      data_[i] = t;
    }
  }

  /// \brief Add element to the end of the vector
  void PushBack(const int& t){
     xsize_ += 1;
     data_.push_back(t);
  }

  /// \brief Return a pointer to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  int* GetArrayPointer() {
    return &(data_[0]);
  }

  /// \brief Return a const point to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  const int* GetConstArrayPointer() const {
    return &(data_[0]);
  }

  // allows access element by element, e.g. this(i) gives data_[i]
  int& operator()(int ix) {return data_[ix];}
  const int& operator()(int ix) const {return data_[ix];}

  /// \brief Insert a given array to the position ix
  /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
  void insert(Array1D<int>& insarr,int ix){
    if (ix<0 || ix>xsize_){
      throw Tantrum("Array1D:insert():: insert index out of bounds.");
    }
    int addsize = insarr.Length();
    xsize_+=addsize;
    int* ptr=insarr.GetArrayPointer();
    data_.insert(data_.begin()+ix,ptr,ptr+addsize);
  }

  /// \brief Insert a given value to the position ix
  /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
  void insert(const int& insval,int ix){
    if (ix<0 || ix>xsize_)
      throw Tantrum("Array1D:insert():: insert index out of bounds.");
    xsize_+=1;
    data_.insert(data_.begin()+ix,insval);
  }

  /// \brief Erase the value from the position ix
  void erase(int ix){
    if (ix<0 || ix>=xsize_)
      throw Tantrum("Array1D:erase():: erase index out of bounds.");
    xsize_-=1;
    data_.erase(data_.begin()+ix);
  }

  /// \brief Dump contents of the array to a file in binary format
  void DumpBinary(FILE* f_out) const {
    fwrite(&xsize_,sizeof(xsize_),1,f_out);
    fwrite(this->GetConstArrayPointer(),sizeof(int),xsize_,f_out);
  }

  /// \brief Read contents of the array from a file in binary format
  void ReadBinary(FILE* f_in){
    fread(&xsize_,sizeof(xsize_),1,f_in);
    data_.resize(xsize_);
    fread(this->GetArrayPointer(),sizeof(int),xsize_,f_in);
  }

  /**********************************************************
  // Methods for interfacing with python 
  **********************************************************/

  // For calling [] in Python
  int& operator[](int i) {return data_[i];}

  /// \brief Dump contents of the array to a file in binary format
  // cannot be read with numpy's fromfile after creation
  void DumpBinary(char *filename){
    FILE *f_out; 
    f_out = fopen(filename,"wb");
    fwrite(&xsize_,sizeof(xsize_),1,f_out);
    fwrite(this->GetConstArrayPointer(),sizeof(int),xsize_,f_out);
    fclose(f_out);
  }

  // read binary file created with DumpBinary
  // Cannot use numpy's from files
  // only for use in c++
  void ReadBinary(char *filename){
    FILE *f_in;
    f_in = fopen(filename,"rb");
    fread(&xsize_,sizeof(xsize_),1,f_in);
    data_.resize(xsize_);
    fread(this->GetArrayPointer(),sizeof(int),xsize_,f_in);
    fclose(f_in);
  }

  // creates binary file that can be read with numpy's fromfile 
  void DumpBinary4py(char *filename){
    ofstream f_out;
    f_out.open(filename, ios::out | ios::binary);
    f_out.write((char*)this->GetArrayPointer(),xsize_*sizeof(int)); // convert array pointer to char string
    f_out.close();
  }

  // can read in DumpBinary4py output, but needs size of vector
  // fromfile can automatically detect size file, so, if need by, one can use numpy's fromfile to determine # of elements
  void ReadBinary4py(char *filename, int n){
    xsize_ = n; 
    ifstream f_in;
    f_in.open(filename, ios::in | ios::binary);
    f_in.read((char*)this->GetArrayPointer(),xsize_*sizeof(int)); // convert array pointer to char string
    f_in.close();
  }

  // Set user-defined list to data_ vector
  // This will work even for string type
  void setArray(vector<int> inarray){
    data_ = inarray; 
    xsize_ = inarray.size();
  }

  // // Sets user-defined 1d numpy array to data_ vector
  void setnpintArray(long* inarray, int n){
    xsize_ = n;
    data_.assign(inarray,inarray+n);
  }
    // This is not to be used for a string type
  void getnpintArray(long* outarray, int n){
    // xsize_ = n;
    // data_.assign(inarray,inarray+n);
    copy(data_.begin(), data_.end(), outarray);
  }

  // Returns data_ vector as a list in python
  // Also acts as a print to see individual elements
  vector<int> flatten(){
    return data_;
  }

  string type(){
    return "int";
  }

};

template<>
class Array1D <double> {
private:
  int xsize_; // private size of vector
public:
  vector<double> data_; // private copy of data vector


  /// Default constructor, which does not allocate any memory
  Array1D(): xsize_(0) {};

  /// \brief Constructor that allocates the memory
  Array1D(const int& nx): xsize_(nx) {
    data_.resize(xsize_);
  }

  /// Constructor that allocates and initializes the data to a value t
  Array1D(const int& nx, const double& t): xsize_(nx) {
    data_.resize(xsize_, t);
  }

  /// \brief Assignment operator copies the data structure by value
  Array1D& operator=(const Array1D &obj) {
    xsize_ = obj.xsize_;
    data_ = obj.data_;
    return *this;
  }

  /// \brief Copy constructor
  Array1D(const Array1D &obj): xsize_(obj.xsize_), data_(obj.data_) {};

  /// Destructor that frees up the memory
  ~Array1D() {data_.clear();}

  /// \brief Function to clear the memory
  void Clear() {
    xsize_ = 0;
    data_.clear();
  }

  /// \brief Returns size in the x-direction
  int XSize() const {return xsize_;}

  /// Returns length (i.e. size in the x-direction)
  int Length() const {return xsize_;}

  /// \brief Resizes the array
  void Resize(const int& nx) {
    xsize_ = nx;
    data_.resize(xsize_);
  }

  /// \brief Resizes the array and sets ALL entries to the specified value
  /// \warning All original data will get lost if this function is used!
  void Resize(const int& nx, const double& t) {
    data_.clear();
    xsize_ = nx;
    data_.resize(xsize_, t);
  }

  /// \brief Set all values in the array to the given value
  void SetValue(const double& t){
    for(int i=0; i < (int)data_.size(); i++){
      data_[i] = t;
    }
  }

  /// \brief Add element to the end of the vector
  void PushBack(const double& t){
     xsize_ += 1;
     data_.push_back(t);
  }

  /// \brief Return a pointer to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  double* GetArrayPointer() {
    return &(data_[0]);
  }

  /// \brief Return a const point to the first element of the data in the
  /// vector so we can use it access the data in array format (e.g. for
  /// passing it to a Fortran program).
  const double* GetConstArrayPointer() const {
    return &(data_[0]);
  }

  // allows access element by element, e.g. this(i) gives data_[i]
  double& operator()(int ix) {return data_[ix];}
  const double& operator()(int ix) const {return data_[ix];}

  /// \brief Insert a given array to the position ix
  /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
  void insert(Array1D<double>& insarr,int ix){
    if (ix<0 || ix>xsize_){
      throw Tantrum("Array1D:insert():: insert index out of bounds.");
    }
    int addsize = insarr.Length();
    xsize_+=addsize;
    double* ptr=insarr.GetArrayPointer();
    data_.insert(data_.begin()+ix,ptr,ptr+addsize);
  }

  /// \brief Insert a given value to the position ix
  /// \note ix=0 means insert at the beginning, ix=xsize_ means insert at the end
  void insert(const double& insval,int ix){
    if (ix<0 || ix>xsize_)
      throw Tantrum("Array1D:insert():: insert index out of bounds.");
    xsize_+=1;
    data_.insert(data_.begin()+ix,insval);
  }

  /// \brief Erase the value from the position ix
  void erase(int ix){
    if (ix<0 || ix>=xsize_)
      throw Tantrum("Array1D:erase():: erase index out of bounds.");
    xsize_-=1;
    data_.erase(data_.begin()+ix);
  }

  /// \brief Dump contents of the array to a file in binary format
  void DumpBinary(FILE* f_out) const {
    fwrite(&xsize_,sizeof(xsize_),1,f_out);
    fwrite(this->GetConstArrayPointer(),sizeof(double),xsize_,f_out);
  }

  /// \brief Read contents of the array from a file in binary format
  void ReadBinary(FILE* f_in){
    fread(&xsize_,sizeof(xsize_),1,f_in);
    data_.resize(xsize_);
    fread(this->GetArrayPointer(),sizeof(double),xsize_,f_in);
  }

  /**********************************************************
  // Methods for interfacing with python 
  **********************************************************/

  // For calling [] in Python
  double& operator[](int i) {return data_[i];}

  /// \brief Dump contents of the array to a file in binary format
  // cannot be read with numpy's fromfile after creation
  void DumpBinary(char *filename){
    FILE *f_out; 
    f_out = fopen(filename,"wb");
    fwrite(&xsize_,sizeof(xsize_),1,f_out);
    fwrite(this->GetConstArrayPointer(),sizeof(double),xsize_,f_out);
    fclose(f_out);
  }

  // read binary file created with DumpBinary
  // Cannot use numpy's from files
  // only for use in c++
  void ReadBinary(char *filename){
    FILE *f_in;
    f_in = fopen(filename,"rb");
    fread(&xsize_,sizeof(xsize_),1,f_in);
    data_.resize(xsize_);
    fread(this->GetArrayPointer(),sizeof(double),xsize_,f_in);
    fclose(f_in);
  }

  // creates binary file that can be read with numpy's fromfile 
  void DumpBinary4py(char *filename){
    ofstream f_out;
    f_out.open(filename, ios::out | ios::binary);
    f_out.write((char*)this->GetArrayPointer(),xsize_*sizeof(double)); // convert array pointer to char string
    f_out.close();
  }

  // can read in DumpBinary4py output, but needs size of vector
  // fromfile can automatically detect size file, so, if need by, one can use numpy's fromfile to determine # of elements
  void ReadBinary4py(char *filename, int n){
    xsize_ = n; 
    ifstream f_in;
    f_in.open(filename, ios::in | ios::binary);
    f_in.read((char*)this->GetArrayPointer(),xsize_*sizeof(double)); // convert array pointer to char string
    f_in.close();
  }

  // Set user-defined list to data_ vector
  // This will work even for string type
  void setArray(vector<double> inarray){
    data_ = inarray; 
    xsize_ = inarray.size();
  }

  // Sets user-defined 1d numpy array to data_ vector
  // This is not to be used for a string type
  void setnpdblArray(double* inarray, int n){
    xsize_ = n;
    data_.assign(inarray,inarray+n);
  }
  // Sets user-defined 1d numpy array to data_ vector
  // This is not to be used for a string type
  void getnpdblArray(double* outarray, int n){
    // xsize_ = n;
    // data_.assign(inarray,inarray+n);
    copy(data_.begin(), data_.end(), outarray);
  }

  // Returns data_ vector as a list in python
  // Also acts as a print to see individual elements
  vector<double> flatten(){
    return data_;
  }

  string type(){
    return "double";
  }
};

#endif /* ARRAY1D_H_SEEN */
