#include "AllClassWrappers.h"

#include "MUQ/Modeling/ModPiece.h"

#include "MUQ/Modeling/Distributions/PyDistribution.h"
#include "MUQ/Modeling/Distributions/Density.h"
#include "MUQ/Modeling/Distributions/DensityProduct.h"
#include "MUQ/Modeling/Distributions/RandomVariable.h"
#include "MUQ/Modeling/Distributions/Distribution.h"
#include "MUQ/Modeling/Distributions/Gaussian.h"
#include "MUQ/Modeling/Distributions/InverseGamma.h"
#include "MUQ/Modeling/Distributions/UniformBox.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Modeling::PythonBindings;
using namespace muq::Modeling;
namespace py = pybind11;

class PyDistributionTramp : public PyDistribution {
public:
  // Inherit the constructors
  using PyDistribution::PyDistribution;

  virtual inline Eigen::VectorXd SampleImpl(std::vector<Eigen::VectorXd> const& inputs) override {
    PYBIND11_OVERLOAD_PURE(Eigen::VectorXd, PyDistribution, SampleImpl, inputs);
  }

  virtual inline double LogDensityImpl(std::vector<Eigen::VectorXd> const& inputs) override {
    PYBIND11_OVERLOAD_PURE(double, PyDistribution, LogDensityImpl, inputs);
  }
private:
};

class PyGaussianTramp : public PyGaussianBase {
public:
    /* Inherit the constructors */
    using PyGaussianBase::PyGaussianBase;

    /* Trampoline (need one for each virtual function) */
    unsigned int Dimension() const override {
      PYBIND11_OVERLOAD(unsigned int,GaussianBase,Dimension);
    }

    Eigen::MatrixXd ApplyCovariance(Eigen::Ref<const Eigen::MatrixXd> const& x) const override{
      PYBIND11_OVERLOAD_PURE(Eigen::MatrixXd, PyGaussianBase, ApplyCovariance, x);
    }

    Eigen::MatrixXd ApplyPrecision(Eigen::Ref<const Eigen::MatrixXd> const& x) const override{
      PYBIND11_OVERLOAD_PURE(Eigen::MatrixXd, PyGaussianBase, ApplyPrecision, x);
    }

    virtual Eigen::MatrixXd ApplyCovSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const override{
      PYBIND11_OVERLOAD_PURE(Eigen::MatrixXd, PyGaussianBase, ApplyCovSqrt, x);
    }

    Eigen::MatrixXd ApplyPrecSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const override{
      PYBIND11_OVERLOAD_PURE(Eigen::MatrixXd, PyGaussianBase, ApplyPrecSqrt, x);
    }

    Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& params) override{
      PYBIND11_OVERLOAD(Eigen::VectorXd, PyGaussianBase, SampleImpl, params);
    }

    void ResetHyperparameters(ref_vector<Eigen::VectorXd> const& params) override{
      PYBIND11_OVERLOAD(void, PyGaussianBase, ResetHyperparameters, params);
    }

    double LogDeterminant() const override{
      PYBIND11_OVERLOAD(double, PyGaussianBase, LogDeterminant);
    }

};

class Publicist : public PyDistribution {
public:
    // Expose protected functions
    using PyDistribution::SampleImpl;
    using PyDistribution::LogDensityImpl;
};

void muq::Modeling::PythonBindings::DistributionWrapper(py::module &m)
{
    py::class_<Distribution, std::shared_ptr<Distribution>> dist(m, "Distribution");
    dist
      .def("LogDensity", (double (Distribution::*)(std::vector<Eigen::VectorXd> const&)) &Distribution::LogDensity)
      .def("LogDensity", (double (Distribution::*)(Eigen::VectorXd const&)) &Distribution::LogDensity)
      .def("GradLogDensity", (Eigen::VectorXd (Distribution::*)(unsigned int, std::vector<Eigen::VectorXd> const&)) &Distribution::GradLogDensity)
      //.def("GradLogDensity", (Eigen::VectorXd (Distribution::*)(unsigned int, Eigen::VectorXd const&)) &Distribution::GradLogDensity)
      .def("Sample", (Eigen::VectorXd (Distribution::*)(std::vector<Eigen::VectorXd> const&)) &Distribution::Sample)
      //.def("Sample", (Eigen::VectorXd (Distribution::*)(Eigen::VectorXd const&)) &Distribution::Sample)
      .def("Sample", (Eigen::VectorXd (Distribution::*)()) &Distribution::Sample)
      .def("AsDensity", &Distribution::AsDensity)
      .def("AsVariable", &Distribution::AsVariable)
      .def_readonly("varSize", &Distribution::varSize)
      .def_readonly("hyperSizes", &Distribution::hyperSizes);


    py::class_<PyDistribution, PyDistributionTramp, Distribution, std::shared_ptr<PyDistribution>> pydist(m, "PyDistribution");
    pydist.def(py::init<unsigned int>());
    pydist.def(py::init<unsigned int, Eigen::VectorXi const&>());
    pydist.def("SampleImpl", (Eigen::VectorXd (PyDistribution::*)(std::vector<Eigen::VectorXd> const&)) &Publicist::SampleImpl);
    pydist.def("LogDensityImpl", (double (PyDistribution::*)(std::vector<Eigen::VectorXd> const&)) &Publicist::LogDensityImpl);

    py::class_<DensityBase, Distribution, ModPiece, std::shared_ptr<DensityBase>> densBase(m, "DensityBase");
    densBase
      .def(py::init<Eigen::VectorXi>());

    py::class_<Density, DensityBase, std::shared_ptr<Density>> dens(m, "Density");
    dens
      .def("GetDistribution", &Density::GetDistribution);

    py::class_<DensityProduct, DensityBase, std::shared_ptr<DensityProduct>> densProd(m,"DensityProduct");
    densProd
      .def(py::init<int>());

    py::class_<RandomVariable, Distribution, ModPiece, std::shared_ptr<RandomVariable>> rv(m, "RandomVariable");
    rv
      .def(py::init<std::shared_ptr<Distribution>>());

    py::class_<UniformBox, Distribution, std::shared_ptr<UniformBox>> uBox(m, "UniformBox");
    uBox
      .def(py::init<Eigen::MatrixXd const&>());

    py::class_<GaussianBase, Distribution, std::shared_ptr<GaussianBase>>(m,"GaussianBase")
      .def("Dimension", &GaussianBase::Dimension)
      .def("GetMean", &GaussianBase::GetMean)
      .def("SetMean", &GaussianBase::SetMean);

    py::class_<Gaussian, GaussianBase, std::shared_ptr<Gaussian>> gauss(m,"Gaussian");
    gauss
      .def(py::init<unsigned int>())
      .def(py::init<unsigned int, Gaussian::InputMask>())
      .def(py::init<Eigen::VectorXd const&>())
      .def(py::init<Eigen::VectorXd const&, Gaussian::InputMask>())
      .def(py::init<Eigen::VectorXd const&, Eigen::MatrixXd const&>())
      .def(py::init<Eigen::VectorXd const&, Eigen::MatrixXd const&, Gaussian::Mode>())
      .def(py::init<Eigen::VectorXd const&, Eigen::MatrixXd const&, Gaussian::Mode, Gaussian::InputMask>())
      .def("GetMode", &Gaussian::GetMode)
      .def("GetCovariance", &Gaussian::GetCovariance)
      .def("GetPrecision", &Gaussian::GetPrecision)
      .def("ApplyPrecision", &Gaussian::ApplyPrecision)
      .def("ApplyCovariance", &Gaussian::ApplyCovariance)
      .def("ApplyCovSqrt", &Gaussian::ApplyCovSqrt)
      .def("ApplyPrecSqrt", &Gaussian::ApplyPrecSqrt)
      .def("SetCovariance", &Gaussian::SetCovariance)
      .def("SetPrecision", &Gaussian::SetPrecision)
      .def("Condition", &Gaussian::Condition);

    py::class_<PyGaussianBase, PyGaussianTramp, GaussianBase, Distribution, std::shared_ptr<PyGaussianBase>>(m, "PyGaussianBase")
        .def(py::init<unsigned int>())
        .def(py::init<unsigned int, Eigen::VectorXi>())
        .def(py::init<Eigen::VectorXd>())
        .def(py::init<Eigen::VectorXd, Eigen::VectorXi>())
        .def("Dimension", &PyGaussianBase::Dimension)
        .def("ApplyCovariance", (Eigen::MatrixXd (PyGaussianBase::*)(Eigen::MatrixXd const&) const) &PyGaussianBase::ApplyCovariance)
        .def("ApplyPrecision", (Eigen::MatrixXd (PyGaussianBase::*)(Eigen::MatrixXd const&) const) &PyGaussianBase::ApplyPrecision)
        .def("ApplyCovSqrt",(Eigen::MatrixXd (PyGaussianBase::*)(Eigen::MatrixXd const&) const)  &PyGaussianBase::ApplyCovSqrt)
        .def("ApplyPrecSqrt", (Eigen::MatrixXd (PyGaussianBase::*)(Eigen::MatrixXd const&) const) &PyGaussianBase::ApplyPrecSqrt)
        .def("LogDeterminant", &PyGaussianBase::LogDeterminant)
        .def("SetMean", &PyGaussianBase::SetMean)
        .def("GetMean", &PyGaussianBase::GetMean)
        .def("SampleImpl", (Eigen::VectorXd (PyGaussianBase::*)(ref_vector<Eigen::VectorXd> const&)) &PyGaussianBase::SampleImpl);


    py::enum_<Gaussian::Mode>(gauss, "Mode")
          .value("Covariance", Gaussian::Mode::Covariance)
          .value("Precision", Gaussian::Mode::Precision)
          .export_values();

    py::enum_<Gaussian::ExtraInputs>(gauss, "ExtraInputs", py::arithmetic())
          .value("None", Gaussian::ExtraInputs::None)
          .value("Mean", Gaussian::ExtraInputs::Mean)
          .value("DiagCovariance", Gaussian::ExtraInputs::DiagCovariance)
          .value("DiagPrecision", Gaussian::ExtraInputs::DiagPrecision)
          .value("FullCovariance", Gaussian::ExtraInputs::FullCovariance)
          .value("FullPrecision", Gaussian::ExtraInputs::FullPrecision)
          .export_values();


    py::class_<InverseGamma, Distribution, std::shared_ptr<InverseGamma>> ig(m, "InverseGamma");
    ig
      .def(py::init<double,double>())
      .def_readonly("alpha", &InverseGamma::alpha)
      .def_readonly("beta", &InverseGamma::beta);
}
