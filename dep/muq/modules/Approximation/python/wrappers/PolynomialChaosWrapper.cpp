#include "AllClassWrappers.h"

#include "MUQ/Approximation/Quadrature/Quadrature.h"
#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"
#include "MUQ/Approximation/PolynomialChaos/PCEFactory.h"
#include "MUQ/Approximation/PolynomialChaos/AdaptiveSmolyakPCE.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/Modeling/ModPiece.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
using namespace muq::Modeling;
using namespace muq::Utilities;
using namespace muq::Approximation;

namespace py = pybind11;

void muq::Approximation::PythonBindings::PolynomialChaosWrapper(py::module &m)
{

  py::class_<PolynomialChaosExpansion, BasisExpansion, std::shared_ptr<PolynomialChaosExpansion>> pce(m,"PolynomialChaosExpansion");
  pce
    .def(py::init<std::shared_ptr<OrthogonalPolynomial> const&, std::shared_ptr<muq::Utilities::MultiIndexSet> const&, Eigen::MatrixXd const&>())
    .def(py::init<std::shared_ptr<OrthogonalPolynomial> const&, std::shared_ptr<muq::Utilities::MultiIndexSet> const&, unsigned int >())
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&, std::shared_ptr<muq::Utilities::MultiIndexSet> const&, Eigen::MatrixXd const&>())
    .def(py::init<std::vector<std::shared_ptr<IndexedScalarBasis>> const&, std::shared_ptr<muq::Utilities::MultiIndexSet> const&, unsigned int>())
    .def("Variance", &PolynomialChaosExpansion::Variance)
    .def("Covariance", &PolynomialChaosExpansion::Covariance)
    .def("Mean", &PolynomialChaosExpansion::Mean)
    .def("Magnitude", &PolynomialChaosExpansion::Magnitude)
    .def("TotalSensitivity", (Eigen::VectorXd (PolynomialChaosExpansion::*)(unsigned int) const) &PolynomialChaosExpansion::TotalSensitivity)
    .def("TotalSensitivity", (Eigen::MatrixXd (PolynomialChaosExpansion::*)() const) &PolynomialChaosExpansion::TotalSensitivity)
    .def("SobolSensitivity", (Eigen::VectorXd (PolynomialChaosExpansion::*)(unsigned int) const) &PolynomialChaosExpansion::SobolSensitivity)
    .def("SobolSensitivity", (Eigen::VectorXd (PolynomialChaosExpansion::*)(std::vector<unsigned int> const&) const) &PolynomialChaosExpansion::SobolSensitivity)
    .def("MainSensitivity", &PolynomialChaosExpansion::MainSensitivity);


  py::class_<PCEFactory, std::shared_ptr<PCEFactory>> pceFactory(m, "PCEFactory");
  pceFactory
    .def(py::init<std::vector<std::shared_ptr<Quadrature>> const&, std::vector<std::shared_ptr<IndexedScalarBasis>> const&>())
    .def(py::init<std::vector<std::shared_ptr<Quadrature>> const&, std::shared_ptr<muq::Utilities::MultiIndex> const&,std::vector<std::shared_ptr<IndexedScalarBasis>> const&>())
    .def(py::init<std::vector<std::shared_ptr<Quadrature>> const&, std::shared_ptr<muq::Utilities::MultiIndex> const&, std::vector<std::shared_ptr<IndexedScalarBasis>> const&,std::shared_ptr<muq::Utilities::MultiIndexSet> const&>())
    .def("Compute", (std::shared_ptr<PolynomialChaosExpansion>(PCEFactory::*)(std::vector<Eigen::VectorXd> const&, std::shared_ptr<MultiIndex> const&) ) &PCEFactory::Compute)
    .def("Compute", (std::shared_ptr<PolynomialChaosExpansion>(PCEFactory::*)(std::shared_ptr<ModPiece> const&) ) &PCEFactory::Compute)
    .def("QuadPts", (std::vector<Eigen::VectorXd>const&(PCEFactory::*)() const) &PCEFactory::QuadPts)
    .def("QuadPts", (std::vector<Eigen::VectorXd>const&(PCEFactory::*)(std::shared_ptr<MultiIndex> const&)) &PCEFactory::QuadPts);


  py::class_<AdaptiveSmolyakPCE, std::shared_ptr<AdaptiveSmolyakPCE>> adaptSmoly(m,"AdaptiveSmolyakPCE");
  adaptSmoly
    .def(py::init<std::shared_ptr<muq::Modeling::ModPiece>, std::vector<std::shared_ptr<Quadrature>>,std::vector<std::shared_ptr<IndexedScalarBasis>>>())
    .def("Compute", [](AdaptiveSmolyakPCE & self,
                       std::shared_ptr<MultiIndexSet> const& fixedSet,
                       py::dict d) {return self.Compute(fixedSet, ConvertDictToPtree(d));})
    .def("Adapt", [](AdaptiveSmolyakPCE & self, py::dict d) {return self.Adapt(ConvertDictToPtree(d));})
    .def("Error", &SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>::Error)
    .def("NumEvals", &SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>::NumEvals)
    .def("ErrorHistory", &SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>::ErrorHistory)
    .def("EvalHistory", &SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>::EvalHistory)
    .def("TimeHistory", &SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>::TimeHistory)
    .def("PointHistory", &SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>::PointHistory)
    .def("TermHistory", &SmolyakEstimator<std::shared_ptr<PolynomialChaosExpansion>>::TermHistory);

}
