#include "AllClassWrappers.h"

#include "MUQ/Approximation/GaussianProcesses/KarhunenLoeveExpansion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

#include "MUQ/Utilities/PyDictConversion.h"

using namespace muq::Utilities;
using namespace muq::Approximation::PythonBindings;
using namespace muq::Approximation;
namespace py = pybind11;

void muq::Approximation::PythonBindings::KLWrapper(py::module &m)
{
  // KarhunenLoeveExpansion class
  py::class_<KarhunenLoeveExpansion, std::shared_ptr<KarhunenLoeveExpansion>>
    KLExpansion(m, "KarhunenLoeveExpansion");
  KLExpansion
    .def(py::init<std::shared_ptr<KernelBase>, Eigen::MatrixXd const&,
                  Eigen::VectorXd const&>())
    .def(py::init( [](std::shared_ptr<KernelBase> const& kern, Eigen::MatrixXd const& pts,
                  Eigen::VectorXd const& wts, py::dict const& d) { return new KarhunenLoeveExpansion(kern, pts, wts, ConvertDictToPtree(d)); } ))
    .def("GetModes", &KarhunenLoeveExpansion::GetModes)
    .def("Evaluate", &KarhunenLoeveExpansion::Evaluate);
}
