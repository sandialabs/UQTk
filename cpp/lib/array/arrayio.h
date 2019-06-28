/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
/// \file arrayio.h
/// \brief Header file for array read/write utilities

#ifndef ARRAYIO_H
#define ARRAYIO_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "Array1D.h"
#include "Array2D.h"


/// \brief Read a datafile from filename in a matrix form
/// and store it in the 2d array data of typename T
/// \note The array data needs to have the correct sizes
template <typename T> void read_datafile(Array2D<T> &data, const char *filename);

/// \brief Read a datafile from filename in a matrix form
/// and store it in the 2d array data if typename T
/// \note The array data is resized to match the file contents
/// \note This function makes two passes: the first pass figures the no. or rows and columns,
/// then the data array is appropriately resized, and the filename is read during second pass
template <typename T> void read_datafileVS(Array2D<T> &data, const char *filename);

/// \brief Read a datafile from filename in a matrix form
/// and store it in a std::vector in column-major storage scheme
/// \note The vector is resized to match the file contents
/// \note This function makes two passes: the first pass figures the no. or rows and columns,
/// then the data vector is appropriately resized, and the filename is read during second pass
template <typename T> void read_datafileVS(std::vector<T> &data, int &nrows, int &ncols, const char *filename);

/// \brief Read a data from filename in a vector form
/// and store it in a 1d array data of typename T
/// \note The array data needs to have the correct size
template <typename T> void read_datafile_1d(Array1D<T>& data, const char* filename);

/// \brief Read a datafile from filename in a vector form
/// and store it in the 1d array data of typename T
/// \note The array data is resized to match the file contents
/// \note This function makes two passes: the first pass figures the no. or rows and columns,
/// then the data array is appropriately resized, and the filename is read during second pass
template <typename T> void read_datafileVS(Array1D<T> &data, const char *filename);

/// \brief Write/append the contents of a 2d array data of typename T
/// to file filename in a matrix form
template <typename T>
void write_datafile(const Array2D<T> &data, const char *filename, const char *action);

/// \brief Write the contents of a 2d array data of typename T
/// to file filename in a matrix form
template <typename T>
void write_datafile(const Array2D<T> &data, const char *filename);

/// \brief Write the contents of a vector of typename T
/// to file filename in a matrix form
template <typename T>
void write_datafile(const std::vector<T> &data, const int &nrows, const int &ncols, const char *storage, const char *filename, const char *action);

/// \brief Write to file filename the number of rows and number of
/// columns on the first line, followed by the contents of a 2d array
/// data of typename T in a matrix form.
template <typename T>
void write_datafile_size(const Array2D<T> &data, const char *filename);

/// \brief Write the contents of a 1d array data of typename T
/// to file filename in a vector form
template <typename T> void write_datafile_1d(const Array1D<T>& data, const char* filename);

#endif // ARRAYIO_H
