#include "AllClassWrappers.h"

#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/MultiLogisticLikelihood.h"
#include "MUQ/Modeling/PyModPiece.h"
#include "MUQ/Modeling/ReplicateOperator.h"
#include "MUQ/Modeling/WorkGraph.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Modeling;
namespace py = pybind11;

class PyModPieceTramp : public PyModPiece {
public:
  /* Inherit the constructors */
  using PyModPiece::PyModPiece;

  /* Trampoline (need one for each virtual function) */
  void EvaluateImpl(std::vector<Eigen::VectorXd> const& input) override {
    PYBIND11_OVERLOAD_PURE(void, PyModPiece, EvaluateImpl, input);
  }

  void GradientImpl(unsigned int                 const  outputDimWrt,
                    unsigned int                 const  inputDimWrt,
                    std::vector<Eigen::VectorXd> const& input,
                    Eigen::VectorXd              const& sensitivity) override {
    PYBIND11_OVERLOAD(void, PyModPiece, GradientImpl, outputDimWrt, inputDimWrt, input, sensitivity);
  }

  void JacobianImpl(unsigned int                 const  outputDimWrt,
                    unsigned int                 const  inputDimWrt,
                    std::vector<Eigen::VectorXd> const& input) override {
    PYBIND11_OVERLOAD(void, PyModPiece, JacobianImpl, outputDimWrt, inputDimWrt, input);
  }

  void ApplyJacobianImpl(unsigned int                 const  outputDimWrt,
                         unsigned int                 const  inputDimWrt,
                         std::vector<Eigen::VectorXd> const& input,
                         Eigen::VectorXd              const& vec) override {
    PYBIND11_OVERLOAD(void, PyModPiece, ApplyJacobianImpl, outputDimWrt, inputDimWrt, input, vec);
  }
};

class Publicist : public PyModPiece {
public:
    // Expose protected functions
    using PyModPiece::EvaluateImpl;
    using PyModPiece::GradientImpl;
    using PyModPiece::JacobianImpl;
    using PyModPiece::ApplyJacobianImpl;

    // Expose protected member variables
    using PyModPiece::outputs;
    using PyModPiece::gradient;
    using PyModPiece::jacobian;
    using PyModPiece::jacobianAction;
};

void muq::Modeling::PythonBindings::ModPieceWrapper(py::module &m)
{
  // Define some functions from the WorkPiece base class
  py::class_<ModPiece, WorkPiece, std::shared_ptr<ModPiece>> mp(m, "ModPiece");
  mp
    .def("Evaluate", (std::vector<Eigen::VectorXd> const& (ModPiece::*)(std::vector<Eigen::VectorXd> const&)) &ModPiece::Evaluate)
    .def("Evaluate", (std::vector<Eigen::VectorXd> const& (ModPiece::*)()) &ModPiece::Evaluate)
    .def_readonly("inputSizes", &ModPiece::inputSizes)
    .def_readonly("outputSizes", &ModPiece::outputSizes)
    .def("GetRunTime", &ModPiece::GetRunTime)
    .def("ResetCallTime", &ModPiece::ResetCallTime)
    .def("GetNumCalls", &ModPiece::GetNumCalls)
    .def("Gradient", (Eigen::VectorXd const& (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::Gradient)
    .def("Jacobian", (Eigen::MatrixXd const& (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&)) &ModPiece::Jacobian)
    .def("ApplyJacobian", (Eigen::VectorXd const& (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::ApplyJacobian)
    .def("GradientByFD", (Eigen::VectorXd (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::GradientByFD)
    .def("JacobianByFD", (Eigen::MatrixXd (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&)) &ModPiece::JacobianByFD)
    .def("ApplyJacobianByFD", (Eigen::VectorXd (ModPiece::*)(unsigned int, unsigned int, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &ModPiece::ApplyJacobianByFD);

  py::class_<PyModPiece, PyModPieceTramp, ModPiece, std::shared_ptr<PyModPiece>> pymp(m, "PyModPiece");
  pymp
    .def(py::init<Eigen::VectorXi const&, Eigen::VectorXi const&>())
    .def("EvaluateImpl", (void (PyModPiece::*)(std::vector<Eigen::VectorXd> const&)) &Publicist::EvaluateImpl)
    .def_readwrite("outputs", &Publicist::outputs)
    .def("GradientImpl", (void (PyModPiece::*)(unsigned int const, unsigned int const, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &Publicist::GradientImpl)
    .def_readwrite("gradient", &Publicist::gradient)
    .def("JacobianImpl", (void (PyModPiece::*)(unsigned int const, unsigned int const, std::vector<Eigen::VectorXd> const&)) &Publicist::JacobianImpl)
    .def_readwrite("jacobian", &Publicist::jacobian)
    .def("ApplyJacobianImpl", (void (PyModPiece::*)(unsigned int const, unsigned int const, std::vector<Eigen::VectorXd> const&, Eigen::VectorXd const&)) &Publicist::ApplyJacobianImpl)
    .def_readwrite("jacobianAction", &Publicist::jacobianAction);

  py::class_<ConstantVector, ModPiece, WorkPiece, std::shared_ptr<ConstantVector>> cv(m, "ConstantVector");
  cv
    .def(py::init<Eigen::VectorXd const&>())
    .def("SetValue",  &ConstantVector::SetValue);

  py::class_<ReplicateOperator, ModPiece, WorkPiece, std::shared_ptr<ReplicateOperator>> ro(m, "ReplicateOperator");
  ro
    .def(py::init<unsigned int, unsigned int>());

  py::class_<ModGraphPiece, ModPiece, WorkPiece, std::shared_ptr<ModGraphPiece>> mgp(m, "ModGraphPiece");
  mgp
    .def("GetGraph", &ModGraphPiece::GetGraph)
    .def("GetConstantPieces", &ModGraphPiece::GetConstantPieces);

  py::class_<MultiLogisticLikelihood, ModPiece, WorkPiece, std::shared_ptr<MultiLogisticLikelihood>> mll(m, "MultiLogisticLikelihood");
  mll
    .def(py::init<unsigned int, Eigen::VectorXi>());
}
