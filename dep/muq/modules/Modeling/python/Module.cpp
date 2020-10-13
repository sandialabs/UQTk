#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::Modeling::PythonBindings;
namespace py = pybind11;

PYBIND11_MODULE(pymuqModeling_, m) {
  WorkPieceWrapper(m);
  ModPieceWrapper(m);
  DistributionWrapper(m);
  CwiseUnaryOperatorsWrapper(m);
  LinearOperatorWrapper(m);
  ODEWrapper(m);
}
