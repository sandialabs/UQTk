//=====================================================================================
//
//                      The UQ Toolkit (UQTk) version 3.1.3
//                          Copyright (2023) NTESS
//                        https://www.sandia.gov/UQToolkit/
//                        https://github.com/sandialabs/UQTk
//
//     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
//     retains certain rights in this software.
//
//     This file is part of The UQ Toolkit (UQTk)
//
//     UQTk is open source software: you can redistribute it and/or modify
//     it under the terms of BSD 3-Clause License
//
//     UQTk is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     BSD 3 Clause License for more details.
//
//     You should have received a copy of the BSD 3 Clause License
//     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
//
//     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
//     Sandia National Laboratories, Livermore, CA, USA
//=====================================================================================
#include <pybind11/pybind11.h>

#include "bcs.h"
#include "Array1D.h"
#include "Array2D.h"

PYBIND11_MODULE(_bcs,m){
  m.def("WBCS",&WBCS);
  //m.def("BCS",static_cast<void (*)(Array2D<double> &, Array1D<double> &, Array1D<double> &, double, Array1D<double> &, int, int, double, int, Array1D<double> &, Array1D<int> &, Array1D<double> &, Array1D<double> &, Array1D<double> &, Array1D<double> &)>(&BCS));
  m.def("BCS",static_cast<void (*)(Array2D<double> &, Array1D<double> &, double &, double, Array1D<double> &, int, int, double, int, Array1D<double> &, Array1D<int> &, Array1D<double> &, Array1D<double> &, Array1D<double> &, double &)>(&BCS));
}
