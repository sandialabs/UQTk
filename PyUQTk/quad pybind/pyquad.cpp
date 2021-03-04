#include <pybind/pybind11.h>

#include "../../cpp/lib/quad/quad.h"

using Quad = Quad::Quad;

namespace py = pybind11;

PYBIND11_MODULE(quad, m) {
    py::class_<Quad>(m,"Quad")
      .def(init<>())
      .def(init<char *,char *,int,int,double,double>())
      .def(init<Array1D<string>&, char *, Array1D<int>&,Array1D<double>&,Array1D<double>&>())
      .def("init",&Quad::init)
      .def("SetAlpha",&Quad:SetAlpha)
      .def("SetBeta",&Quad::SetBeta)
      .def("SetDomain",[](Quad quad, Array1D<double>& a, Array1D<double>& b){return quad.SetDomain(a,b);})
      .def("SetDomain",[](Quad quad, Array1D<double>& a){return quad.SetDomain(a);})
      .def("GetDomain",[](Quad quad, Array1D<double>& a, Array1D<double>& b){return quad.GetDomain(a,b);})
      .def("GetDomain",[](Quad quad, Array1D<double>& a){return quad.GetDomain(a);})
      .def("SetRule",[](Quad quad,Array2D<double>& q, Array1D<double>& w){return quad.SetRule(q,w)})
      .def("SetRule",[](Quad quad,Array2D<double>& q, Array1D<double>& w, Array2D<int>& ind){return quad.SetRule(q,w,ind)})
      .def("SetRule",[](Quad quad){return quad.SetRule()})
      .def("SetQdpts",&Quad::SetQdpts)
      .def("SetWghts",&Quad::SetWghts)
      .def("SetLevel",&Quad::SetLevel)
      .def("nextLevel",&Quad::nextLevel)
      .def("GetNQ",&Quad::GetNQ)
      .def("SetVerbosity",&Quad::SetVerbosity)
      .def_property_readonly("aa_",&Quad::aa_)
      .def_property_readonly("bb_",&Quad::aa_)

}
