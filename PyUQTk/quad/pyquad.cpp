#include <pybind11/pybind11.h>

#include "quad.h"

namespace py = pybind11;

PYBIND11_MODULE(_quad, m) {
    py::class_<Quad>(m,"Quad")
      .def(pybind11::init<>())
      .def(pybind11::init<char *,char *,int,int,double,double>())
      .def(pybind11::init<Array1D<string>&, char *, Array1D<int>&,Array1D<double>&,Array1D<double>&>())
      .def("init",&Quad::init)
      .def("SetAlpha",&Quad::SetAlpha)
      .def("SetBeta",&Quad::SetBeta)
      .def("SetDomain",static_cast<void (Quad::*)(Array1D<double>&,Array1D<double>&)>(&Quad::SetDomain))
      .def("SetDomain",static_cast<void (Quad::*)(Array1D<double>&)>(&Quad::SetDomain))
      .def("GetDomain",static_cast<void (Quad::*)(Array1D<double>&,Array1D<double>&) const>(&Quad::GetDomain))
      .def("GetDomain",static_cast<void (Quad::*)(Array1D<double>&) const>(&Quad::GetDomain))
      .def("SetRule",static_cast<void (Quad::*)(Array2D<double>&,Array1D<double>&)>(&Quad::SetRule))
      .def("SetRule",static_cast<void (Quad::*)(Array2D<double>&,Array1D<double>&,Array2D<int>&)>(&Quad::SetRule))
      .def("SetRule",static_cast<void (Quad::*)()>(&Quad::SetRule))
      .def("SetQdpts",&Quad::SetQdpts)
      .def("SetWghts",&Quad::SetWghts)
      .def("SetLevel",&Quad::SetLevel)
      .def("nextLevel",&Quad::nextLevel)
      .def("GetNQ",&Quad::GetNQ)
      .def("SetVerbosity",&Quad::SetVerbosity)
      ;
}
