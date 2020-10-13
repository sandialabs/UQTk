#include "AllClassWrappers.h"

#include "MUQ/Approximation/Polynomials/BasisExpansion.h"
#include "MUQ/Approximation/Polynomials/HermiteFunction.h"
#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"
#include "MUQ/Approximation/Polynomials/Jacobi.h"
#include "MUQ/Approximation/Polynomials/Laguerre.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"
#include "MUQ/Approximation/Polynomials/Monomial.h"
#include "MUQ/Approximation/Polynomials/MonotoneExpansion.h"
#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"
#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"
#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
using namespace muq::Modeling;

namespace py = pybind11;

void muq::Approximation::PythonBindings::PolynomialsWrapper(py::module &m)
{
  py::class_<BasisExpansion, ModPiece, std::shared_ptr<BasisExpansion>> baseExp(m, "BasisExpansion");
  baseExp
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&>())
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&, bool>())
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&, std::shared_ptr<muq::Utilities::MultiIndexSet>>())
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&, std::shared_ptr<muq::Utilities::MultiIndexSet>, bool>())
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&, std::shared_ptr<muq::Utilities::MultiIndexSet>, Eigen::MatrixXd const&>())
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&, std::shared_ptr<muq::Utilities::MultiIndexSet>, Eigen::MatrixXd const&, bool>())
    .def("NumTerms", &BasisExpansion::NumTerms)
    .def("SecondDerivative", (Eigen::MatrixXd (BasisExpansion::*)(unsigned, unsigned, unsigned, Eigen::VectorXd const&, Eigen::MatrixXd const&)) &BasisExpansion::SecondDerivative)
    .def("SecondDerivative", (Eigen::MatrixXd (BasisExpansion::*)(unsigned, unsigned, unsigned, Eigen::VectorXd const&)) &BasisExpansion::SecondDerivative)
    .def("GetCoeffs", &BasisExpansion::GetCoeffs)
    .def("SetCoeffs", &BasisExpansion::SetCoeffs)
    .def("Multis", &BasisExpansion::Multis)
    .def("BuildVandermonde", &BasisExpansion::BuildVandermonde)
    .def("BuildDerivMatrix", &BasisExpansion::BuildDerivMatrix);

  py::class_<IndexedScalarBasis, std::shared_ptr<IndexedScalarBasis>> iSB(m, "IndexedScalarBasis");
  iSB
    .def_static("Construct", &IndexedScalarBasis::Construct)
    .def("EvaluateAllTerms", &HermiteFunction::EvaluateAllTerms);

  py::class_<OrthogonalPolynomial, IndexedScalarBasis, std::shared_ptr<OrthogonalPolynomial>> orthPoly(m, "OrthogonalPolynomial");
  orthPoly
    .def_static("Construct", &OrthogonalPolynomial::Construct)
    .def("BasisEvaluate", &OrthogonalPolynomial::BasisEvaluate)
    .def("EvaluateAllTerms", &OrthogonalPolynomial::DerivativeEvaluate)
    .def("Normalization", &OrthogonalPolynomial::Normalization);

  py::class_<HermiteFunction, IndexedScalarBasis, std::shared_ptr<HermiteFunction>> hermFunc(m, "HermiteFunction");
  hermFunc
    .def(py::init<>())
    .def("BasisEvaluate", &HermiteFunction::BasisEvaluate)
    .def("DerivativeEvaluate", &HermiteFunction::DerivativeEvaluate);

  py::class_<Jacobi, OrthogonalPolynomial, std::shared_ptr<Jacobi>> jacob(m, "Jacobi");
  jacob
    .def(py::init<>())
    .def(py::init<const double>())
    .def(py::init<const double, const double>())
    .def("DerivativeEvaluate", &Jacobi::DerivativeEvaluate)
    .def("Normalization", &Jacobi::Normalization);

  py::class_<Laguerre, OrthogonalPolynomial, std::shared_ptr<Laguerre>> laguerre(m, "Laguerre");
  laguerre
    .def(py::init<>())
    .def(py::init<const double>())
    .def("DerivativeEvaluate", &Laguerre::DerivativeEvaluate)
    .def("Normalization", &Laguerre::Normalization);

  py::class_<Legendre, OrthogonalPolynomial, std::shared_ptr<Legendre>> legend(m, "Legendre");
  legend
    .def(py::init<>())
    .def("DerivativeEvaluate", &Legendre::DerivativeEvaluate)
    .def("Normalization", &Legendre::Normalization);

  py::class_<Monomial, IndexedScalarBasis, std::shared_ptr<Monomial>> mono(m, "Monomial");
  mono
    .def(py::init<>())
    .def("BasisEvaluate", &Monomial::BasisEvaluate)
    .def("DerivativeEvaluate", &Monomial::DerivativeEvaluate);

  py::class_<MonotoneExpansion, ModPiece, std::shared_ptr<MonotoneExpansion>> monotone(m, "MonotoneExpansion");
  monotone
    .def(py::init<std::shared_ptr<BasisExpansion>>())
    .def(py::init<std::shared_ptr<BasisExpansion>, bool>())
    .def(py::init<std::vector<std::shared_ptr<BasisExpansion>> const&, std::vector<std::shared_ptr<BasisExpansion>> const&>())
    .def(py::init<std::vector<std::shared_ptr<BasisExpansion>> const&, std::vector<std::shared_ptr<BasisExpansion>> const&, bool>())
    .def("NumTerms", &MonotoneExpansion::NumTerms)
    .def("Head", &MonotoneExpansion::Head)
    .def("EvaluateInverse", (Eigen::VectorXd (MonotoneExpansion::*)(Eigen::VectorXd const&) const) &MonotoneExpansion::EvaluateInverse)
    .def("EvaluateInverse", (Eigen::VectorXd (MonotoneExpansion::*)(Eigen::VectorXd const&, Eigen::VectorXd const&) const) &MonotoneExpansion::EvaluateInverse)
    .def("EvaluateForward", &MonotoneExpansion::EvaluateForward)
    .def("GetCoeffs", &MonotoneExpansion::GetCoeffs)
    .def("SetCoeffs", &MonotoneExpansion::SetCoeffs)
    .def("LogDeterminant", (double (MonotoneExpansion::*)(Eigen::VectorXd const&)) &MonotoneExpansion::LogDeterminant)
    .def("LogDeterminant", (double (MonotoneExpansion::*)(Eigen::VectorXd const&, Eigen::VectorXd const&)) &MonotoneExpansion::LogDeterminant)
    .def("GradLogDeterminant", (Eigen::VectorXd (MonotoneExpansion::*)(Eigen::VectorXd const&)) &MonotoneExpansion::GradLogDeterminant)
    .def("GradLogDeterminant", (Eigen::VectorXd (MonotoneExpansion::*)(Eigen::VectorXd const&, Eigen::VectorXd const&)) &MonotoneExpansion::GradLogDeterminant);

  py::class_<PhysicistHermite, OrthogonalPolynomial, std::shared_ptr<PhysicistHermite>> physHerm(m, "PhysicistHermite");
  physHerm
    .def(py::init<>())
    .def("DerivativeEvaluate", &PhysicistHermite::DerivativeEvaluate)
    .def("Normalization", &PhysicistHermite::Normalization);

  py::class_<ProbabilistHermite, OrthogonalPolynomial, std::shared_ptr<ProbabilistHermite>> probHerm(m, "ProbabilistHermite");
  probHerm
    .def(py::init<>())
    .def("DerivativeEvaluate", &ProbabilistHermite::DerivativeEvaluate)
    .def("Normalization", &ProbabilistHermite::Normalization);
}
