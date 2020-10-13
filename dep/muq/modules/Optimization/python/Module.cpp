#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::Optimization::PythonBindings;
namespace py = pybind11;

PYBIND11_MODULE(pymuqOptimization, m) {
  CostFunctionWrapper(m);
  OptimizationWrapper(m);
}
