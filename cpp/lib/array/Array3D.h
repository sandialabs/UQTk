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

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/// \file Array3D.h
/// \brief 3D Array class for any type T


#ifndef ARRAY3D_H_SEEN
#define ARRAY3D_H_SEEN

#include <stddef.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>

using namespace std;

/// \class  Array3D
/// \brief  Stores data of any type T in a 3D array
///
/// This class also provides a Fortran-like access operator ()
/// as well as a function to access the data in the array through a pointer that
/// can be passed to F77 or C routines.
/// \author  Bert Debusschere <bjdebus@sandia.gov>
/// \date  Jan 2005
/// \note  Inspired by Helgi Adalsteinsson's Array class implementation
/// \todo  Define copy constructor
/// \todo Several functions, e.g. insert/erase columns/rows, available in Array1D and Array2D, are missing.
template <typename T>
class Array3D {
  public:
    /// \brief Default constructor, which does not allocate any memory
    Array3D(): xsize_(0), ysize_(0), zsize_(0) {};

    /// \brief Constructor that allocates the memory
    Array3D(const size_t& nx, const size_t& ny, const size_t& nz):
      xsize_(nx), ysize_(ny), zsize_(nz) {
      data_.resize(xsize_*ysize_*zsize_);
    }

    /// \brief Constructor that allocates and initializes the data
    Array3D(const size_t& nx, const size_t& ny, const size_t& nz, const T& t):
      xsize_(nx), ysize_(ny), zsize_(nz) {
      data_.resize(xsize_*ysize_*zsize_ , t);
    }

      /// \brief Destructor that frees up the memory
    ~Array3D() {data_.clear();}

    /// \brief Function to clear the memory
    void Clear() {
      xsize_ = 0;
      ysize_ = 0;
      zsize_ = 0;
      data_.clear();
    }

    /// \brief Returns size in the x-direction
    size_t XSize() const {return xsize_;}
    /// \brief Returns size in the y-direction
    size_t YSize() const {return ysize_;}
    /// \brief Returns size in the z-direction
    size_t ZSize() const {return zsize_;}

    /// \brief Resizes the array
    /// \warning In its current implementation, most of the original data
    /// will get lost if the xsize or ysize changes as this changes the indexing for all entries.
    /// \todo Write a better implementation that preserves the original data by
    /// copying it to a temporary array and putting the elements back where they were before.
    /// This would bring this resize() command more closely in line with vector::resize()
    /// function in the original vector class.
    void Resize(const size_t& nx, const size_t& ny, const size_t& nz) {
      xsize_ = nx;
      ysize_ = ny;
      zsize_ = nz;
      data_.resize(xsize_*ysize_*zsize_);
    }

    /// \brief Resizes the array and sets ALL entries to the specified value
    /// \warning All original data will get lost if this function is used!
    /// \todo Write an implementation that is more closely follows the resize
    /// command in the vector class, which keeps the original elements and only
    /// initializes the new elements.
    void Resize(const size_t& nx, const size_t& ny, const size_t& nz, const T& t) {
      data_.clear();
      xsize_ = nx;
      ysize_ = ny;
      zsize_ = nz;
      data_.resize(xsize_*ysize_*zsize_ , t);
    }

    /// \brief Set all values in the array to the given value
    void SetValue(const T& t){
      for(size_t i=0; i < data_.size(); i++){
        data_[i] = t;
      }
    }

    /// \brief Return a pointer to the first element of the data in the
    /// vector so we can use it access the data in array format (e.g. for
    /// passing it to a Fortran program).
    T* GetArrayPointer() {
      return &(data_[0]);
    }

    /// \brief Return a const pointer to the first element of the data in the
    /// vector so we can use it access the data in array format (e.g. for
    /// passing it to a Fortran program).
    const T* GetConstArrayPointer() const {
      return &(data_[0]);
    }

    /// \brief Fortran-like () operator to access values in the 3D data array
    ///
    /// If "my_data" is an object of type Array3D, then its array values can
    /// be accessed as my_data(ix,iy,iz), where ix, iy, iz are the indices in the
    /// x, y, and z dimensions respectively.
    T& operator()(size_t ix, size_t iy, size_t iz) {return data_[ix+xsize_*(iy+ysize_*iz)];}

    /// \brief Fortran-like () const operator to access values in the 3D data array
    ///
    /// If "my_data" is an object of type Array3D, then its array values can
    /// be accessed as my_data(ix,iy,iz), where ix, iy, iz are the indices in the
    /// x, y, and z dimensions respectively.
    const T& operator()(size_t ix, size_t iy, size_t iz) const {return data_[ix+xsize_*(iy+ysize_*iz)];}

    /// \brief Dump contents of the array to a file in binary format
    void DumpBinary(FILE* f_out) const {
      fwrite(&xsize_,sizeof(xsize_),1,f_out);
      fwrite(&ysize_,sizeof(ysize_),1,f_out);
      fwrite(&zsize_,sizeof(zsize_),1,f_out);
      fwrite(this->GetConstArrayPointer(),sizeof(T),xsize_*ysize_*zsize_,f_out);
    }
    
    /// \brief Dump contents of the array to a file in text format
    /// Added by Maher Salloum
    /// When post-processing (in matlab for example), one has to transpose each 2-D
    /// sub-matrix imported from the text file.
    void DumpText(std::ofstream& f_out) const {
      vector<double>::const_iterator it1;
      vector<double>::const_iterator it2;
      it2=data_.begin();
      
      for (int iz=0;iz<zsize_;iz++) {
        for (int iy=0;iy<ysize_;iy++) {
          it1=it2;
          advance(it2,xsize_);
          std::copy(it1,it2,std::ostream_iterator<T>(f_out," "));
          f_out << endl;
        }
      }
      
    }

    /// \brief Read contents of the array from a file in binary format
    void ReadText(FILE* f_in){
      fread(&xsize_,sizeof(xsize_),1,f_in);
      fread(&ysize_,sizeof(ysize_),1,f_in);
      fread(&zsize_,sizeof(zsize_),1,f_in);
      data_.resize(xsize_*ysize_*zsize_);
      fread(this->GetArrayPointer(),sizeof(T),xsize_*ysize_*zsize_,f_in);
    }

    /// \brief Read contents of the array from a file in text format
    /// Added by Maher Salloum
    void ReadBinary(std::ifstream& f_in){
      
      typedef std::istream_iterator<T> istream_iterator;
      std::copy(istream_iterator(f_in),istream_iterator(),data_.begin());
    }
       

  private:

    /// \brief Copy constructor, which is made private so it would not be used inadvertently
    /// (until we define a proper copy constructor)
    Array3D(const Array3D &obj) {};

    /// \brief Number of elements in the x-dimension
    size_t xsize_;
    /// \brief Number of elements in the y-dimension
    size_t ysize_;
    /// \brief Number of elements in the z-dimension
    size_t zsize_;

    /// \brief Data in the array with size = xsize_ * ysize_ * zsize_
    ///
    /// The data is stored with the fastest running index in the x-dimension
    /// then the y-dimension and the slowest one in the z-dimension. The indices
    /// in every dimension run from 0 to their respective "size-1"
    vector<T> data_;
};

#endif /* ARRAY3D_H_SEEN */
