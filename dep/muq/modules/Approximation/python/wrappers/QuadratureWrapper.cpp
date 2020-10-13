#include "AllClassWrappers.h"

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Quadrature/ClenshawCurtisQuadrature.h"
#include "MUQ/Approximation/Quadrature/GaussPattersonQuadrature.h"
#include "MUQ/Approximation/Quadrature/ExponentialGrowthQuadrature.h"
#include "MUQ/Approximation/Quadrature/FullTensorQuadrature.h"
#include "MUQ/Approximation/Quadrature/SmolyakQuadrature.h"
#include "MUQ/Approximation/Quadrature/AdaptiveSmolyakQuadrature.h"

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

void muq::Approximation::PythonBindings::QuadratureWrapper(py::module &m)
{

  py::class_<Quadrature, std::shared_ptr<Quadrature>> quadBase(m, "Quadrature");
  quadBase
    .def("Compute",(void (Quadrature::*)(unsigned int)) &Quadrature::Compute)
    .def("Compute",(void (Quadrature::*)(Eigen::RowVectorXi const&)) &Quadrature::Compute)
    .def("Dim", &Quadrature::Dim)
    .def("Points", &Quadrature::Points)
    .def("Weights", &Quadrature::Weights);

  py::class_<GaussQuadrature, Quadrature, std::shared_ptr<GaussQuadrature>> quadFunc(m, "GaussQuadrature");
  quadFunc
    .def(py::init<std::shared_ptr<OrthogonalPolynomial>>())
    .def(py::init<std::shared_ptr<OrthogonalPolynomial>, int>())
    .def("Compute", &GaussQuadrature::Compute);

  py::class_<ClenshawCurtisQuadrature, Quadrature, std::shared_ptr<ClenshawCurtisQuadrature>> ccQuad(m,"ClenshawCurtisQuadrature");
  ccQuad
    .def(py::init<>())
    .def(py::init<bool>())
    .def("Compute", &ClenshawCurtisQuadrature::Compute);

  py::class_<GaussPattersonQuadrature, Quadrature, std::shared_ptr<GaussPattersonQuadrature>> gpQuad(m,"GaussPattersonQuadrature");
  gpQuad
    .def(py::init<>())
    .def("Compute", &GaussPattersonQuadrature::Compute);

  py::class_<ExponentialGrowthQuadrature, Quadrature, std::shared_ptr<ExponentialGrowthQuadrature>> egQuad(m,"ExponentialGrowthQuadrature");
  egQuad
    .def(py::init<std::shared_ptr<Quadrature>>())
    .def("Compute", &ExponentialGrowthQuadrature::Compute);

  py::class_<FullTensorQuadrature, Quadrature, std::shared_ptr<FullTensorQuadrature>> tensQuad(m,"FullTensorQuadrature");
  tensQuad
    .def(py::init<unsigned int, std::shared_ptr<Quadrature>>())
    .def(py::init<unsigned int, std::shared_ptr<Quadrature>, unsigned int>())
    .def(py::init<std::vector<std::shared_ptr<Quadrature>>, Eigen::RowVectorXi>());

  py::class_<SmolyakQuadrature, Quadrature, std::shared_ptr<SmolyakQuadrature>> smolyQuad(m,"SmolyakQuadrature");
  smolyQuad
    .def(py::init<unsigned int, std::shared_ptr<Quadrature> const&>())
    .def(py::init<std::vector<std::shared_ptr<Quadrature>> const&>())
    .def("Compute", (void (SmolyakQuadrature::*)(unsigned int)) &SmolyakQuadrature::Compute)
    .def("Compute", (void (SmolyakQuadrature::*)(Eigen::RowVectorXi const&)) &SmolyakQuadrature::Compute)
    .def("Compute", (void (SmolyakQuadrature::*)(std::shared_ptr<muq::Utilities::MultiIndexSet> const&)) &SmolyakQuadrature::Compute)
    .def("ComputeWeights", &SmolyakQuadrature::ComputeWeights)
    .def("BuildMultis", &SmolyakQuadrature::BuildMultis);

  // py::class_<SmolyakEstimator<Eigen::VectorXd>, std::shared_ptr<SmolyakEstimator<Eigen::VectorXd>>> smolyBase(m,"SmolyakEstimatorVector");
  // smolyBase
  //   .def(py::init<std::shared_ptr<muq::Modeling::ModPiece>>())
  //   .def("Compute", [](SmolyakEstimator<Eigen::VectorXd> & self,
  //                      std::shared_ptr<MultiIndexSet> const& fixedSet,
  //                      py::dict d) {self.Compute(fixedSet, ConvertDictToPtree(d));})
  //   .def("Adapt", [](SmolyakEstimator<Eigen::VectorXd> & self,
  //                   py::dict d) {self.Adapt(ConvertDictToPtree(d));})
  //   .def("Error", &SmolyakEstimator<Eigen::VectorXd>::Error)
  //   .def("NumEvals", &SmolyakEstimator<Eigen::VectorXd>::NumEvals);
  //

  py::class_<AdaptiveSmolyakQuadrature, std::shared_ptr<AdaptiveSmolyakQuadrature>> adaptSmoly(m,"AdaptiveSmolyakQuadrature");
  adaptSmoly
    .def(py::init<std::shared_ptr<muq::Modeling::ModPiece>, std::vector<std::shared_ptr<Quadrature>>>())
    .def("Compute", [](AdaptiveSmolyakQuadrature & self,
                       std::shared_ptr<MultiIndexSet> const& fixedSet,
                       py::dict d) {return self.Compute(fixedSet, ConvertDictToPtree(d));})
    .def("Adapt", [](AdaptiveSmolyakQuadrature & self, py::dict d) {return self.Adapt(ConvertDictToPtree(d));})
    .def("Error", &SmolyakEstimator<Eigen::VectorXd>::Error)
    .def("NumEvals", &SmolyakEstimator<Eigen::VectorXd>::NumEvals);

}
