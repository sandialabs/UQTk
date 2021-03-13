#include <pybind11/pybind11.h>

#include "PCBasis.h"
#include "PCSet.h"

namespace py = pybind11;

PYBIND11_MODULE(_pce, m) {
    py::class_<PCBasis>(m,"PCBasis")
      .def(py::init<const string, const double, const double, const int>())
      .def("Init1dQuadPoints",&PCBasis::Init1dQuadPoints)
      .def("Eval1dBasisAtQuadPoints",&PCBasis::Eval1dBasisAtQuadPoints)
      .def("Eval1dBasisAtCustPoints",&PCBasis::Eval1dBasisAtCustPoints)
      .def("EvalBasis",static_cast<double (PCBasis::*)(const double &, Array1D<double> &) const>(&PCBasis::EvalBasis))
      .def("EvalBasis",static_cast<double (PCBasis::*)(const double &, const int, double *) const>(&PCBasis::EvalBasis))
      .def("Eval1dNormSq_Exact",&PCBasis::Eval1dNormSq_Exact)
      .def("EvalDerivBasis",&PCBasis::EvalDerivBasis)
      .def("Eval1dDerivBasisAtCustPoints",&PCBasis::Eval1dDerivBasisAtCustPoints)
      .def("Eval2ndDerivBasis",&PCBasis::Eval2ndDerivBasis)
      .def("Eval2ndDerivCustPoints",&PCBasis::Eval2ndDerivCustPoints)
      .def("Get1dNormsSq",&PCBasis::Get1dNormsSq)
      .def("Get1dNormsSqExact",&PCBasis::Get1dNormsSqExact)
      .def("GetRandSample",static_cast<void (PCBasis::*)(Array1D<double>&)>(&PCBasis::GetRandSample))
      .def("GetRandSample",static_cast<void (PCBasis::*)(double*, const int&)>(&PCBasis::GetRandSample))
      .def("GetSeed",&PCBasis::GetSeed)
      .def("SeedRandNumGen",&PCBasis::SeedRandNumGen)
      .def("GetQuadRule",&PCBasis::GetQuadRule)
      .def("GetQuadPoints",&PCBasis::GetQuadPoints)
      .def("GetQuadWeights",&PCBasis::GetQuadWeights)
      .def("GetQuadIndices",&PCBasis::GetQuadIndices)
      .def("GetBasisAtQuadPoints",&PCBasis::GetBasisAtQuadPoints)
      .def("GetPCType",&PCBasis::GetPCType)
      .def("GetAlpha",&PCBasis::GetAlpha)
      .def("GetBeta",&PCBasis::GetBeta)
      ;
}
