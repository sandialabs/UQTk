#include "AllClassWrappers.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <string>

#include <functional>
#include <vector>

#include "MUQ/Optimization/NLoptOptimizer.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/Modeling/Python/PyAny.h"

using namespace muq::Utilities;
using namespace muq::Optimization;
namespace py = pybind11;

void PythonBindings::OptimizationWrapper(pybind11::module &m) {
  py::class_<NLoptOptimizer, std::shared_ptr<NLoptOptimizer> > opt(m, "NLoptOptimizer");
  opt.def(py::init( [](std::shared_ptr<CostFunction> cost, py::dict d) { return new NLoptOptimizer(cost, ConvertDictToPtree(d)); }));
  //opt.def("AddInequalityConstraint", &NLoptOptimizer::AddInequalityConstraint);
  opt.def("Solve", (std::pair<Eigen::VectorXd, double> (NLoptOptimizer::*)(std::vector<Eigen::VectorXd> const&)) &NLoptOptimizer::Solve);
}
