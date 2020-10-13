#include "AllClassWrappers.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <string>

#include <functional>
#include <vector>

#include "MUQ/Optimization/ModPieceCostFunction.h"

using namespace muq::Modeling;
using namespace muq::Optimization;
namespace py = pybind11;

void PythonBindings::CostFunctionWrapper(py::module &m) {
  py::class_<CostFunction, std::shared_ptr<CostFunction> > cost(m, "CostFunction");
  cost
    .def("Cost", (double (CostFunction::*)(std::vector<Eigen::VectorXd> const&)) &CostFunction::Cost)
    .def("Gradient", (Eigen::VectorXd const& (CostFunction::*)(unsigned int const, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &CostFunction::Gradient);

  py::class_<ModPieceCostFunction, CostFunction, std::shared_ptr<ModPieceCostFunction> > modCost(m, "ModPieceCostFunction");
  modCost
    .def(py::init( [](std::shared_ptr<ModPiece> cost) { return new ModPieceCostFunction(cost); }));
}
