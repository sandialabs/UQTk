//=====================================================================================
//
//                      The UQ Toolkit (UQTk) version @UQTKVERSION@
//                          Copyright (@UQTKYEAR@) NTESS
//                        https://www.sandia.gov/UQToolkit/
//                        https://github.com/sandialabs/UQTk
//
//     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

#include "kle.h"
#include <iostream>
#include "ftndefs.h"
#include "Array1D.h"
#include "Array2D.h"

namespace py = pybind11;

PYBIND11_MODULE(_kle, m) {
  py::class_<KLDecompUni>(m,"KLDecompUni"){
    .def(py::init<>())
    .def(py::init<const Array1D<double>&>)
    .def("Init",&KLDecompUni::Init)
    .def("SetWeights",static_cast<void (KLDecompUni::*)(const Array1D<double>&)>(&KLDecompUni::SetWeights))
    .def("SetWeights",static_cast<void (KLDecompUni::*)(const double *, const int)>(&KLDecompUni::SetWeights))
    .def("decompose",static_cast<int (KLDecompUni::*)(const Array2D<double>&, const int&)>(&KLDecompUni::decompose))
    .def("decompose",static_cast<int (KLDecompUni::*)(const double*, const int&)>(&KLDecompUni::decompose))
    .def("KLproject",&KLDecompUni::KLproject)
    .def("eigenvalues",static_cast<const Array1D<double> (KLDecompUni::*)() const>(&KLDecompUni::eigenvalues))
    .def("eigenvalues",static_cast<void (KLDecompUni::*)(const int,double *) const>(&KLDecompUni::eigenvalues))
    .def("KLModes",static_cast<const Array2D<double> (KLDecompUni::*)() const>(&KLDecompUni::KLmodes))
    .def("KLModes",static_cast<void (KLDecompUni::*)(const int, const int, double *) const>(&KLDecompUni::eigenvalues))
    .def("meanRealiz",&KLDecompUni::meanRealiz)
    .def("truncRealiz",&KLDecompUni::truncRealiz)
  }

}
