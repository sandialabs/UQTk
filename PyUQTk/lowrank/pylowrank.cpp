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

#include "lowrank.h"
#include "Array1D.h"
#include "Array2D.h"

namespace py = pybind11;

PYBIND11_MODULE(_lowrank, m) {
    py::class_<CanonicalTensor>(m,"CanonicalTensor")
      .def(pybind11::init<>())
      .def(pybind11::init<Array1D<double>&, Array1D<Array2D<double>> >())
      .def("create", static_cast<void (CanonicalTensor::*)(string, int, Array1D<int>&)>(&CanonicalTensor::create))
      .def("getSpaceVec", static_cast< Array1D<double> (CanonicalTensor::*)(int, int)>(&CanonicalTensor::getSpaceVec))
      .def("setSpaceVec", static_cast<void (CanonicalTensor::*)(int, int, Array1D<double>&)>(&CanonicalTensor::setSpaceVec))
      .def("display",static_cast<void (CanonicalTensor::*)()>(&CanonicalTensor::display))
      .def("prodscal", static_cast<double (*)(CanonicalTensor&,CanonicalTensor&)>(&CanonicalTensor::prodscal))
      .def("norm",static_cast<double (CanonicalTensor::*)()>(&CanonicalTensor::norm))
      ;

    py::class_<PCBases>(m,"PCBases")
      .def(pybind11::init<Array1D<string>&, Array1D<int> >())
      .def("functionEval",static_cast<void (PCBases::*)(Array2D<double>&, int, Array2D<double>&)>(&PCBases::functionEval))
      ;

    py::class_<PLBases>(m,"PLBases")
      .def("functionEval", static_cast<void (PLBases::*)(Array2D<double>&, int, Array2D<double>&)>(&PLBases::functionEval))
      ;

    py::class_<FunctionalTensor>(m, "FunctionalTensor")
      .def(pybind11::init<>())
      .def(pybind11::init<CanonicalTensor&>())
      .def(pybind11::init<CanonicalTensor&, Array1D<FunctionalBases*>&>())
      .def("tensorEval", static_cast<Array1D<double> (FunctionalTensor::*)(Array2D<double>&)>(&FunctionalTensor::tensorEval))
      .def("tensorEval", static_cast<Array1D<double> (FunctionalTensor::*)(Array1D<Array2D<double>>&, int)>(&FunctionalTensor::tensorEval))
      ;

    py::class_<CanonicalTensorALSLeastSquares>(m, "CanonicalTensorALSLeastSquares")
      .def(pybind11::init<>())
      .def("solveDirect", static_cast<CanonicalTensor (CanonicalTensorALSLeastSquares::*)(Array2D<double>&, Array1D<double>&, bool)>(&CanonicalTensorALSLeastSquares::solveDirect))
      ;
}
