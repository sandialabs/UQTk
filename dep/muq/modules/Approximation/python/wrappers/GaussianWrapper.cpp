#include "AllClassWrappers.h"

#include "MUQ/Approximation/GaussianProcesses/GaussianProcess.h"
#include "MUQ/Approximation/GaussianProcesses/ObservationInformation.h"
#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
using namespace muq::Approximation;
namespace py = pybind11;

void muq::Approximation::PythonBindings::GaussianWrapper(py::module &m)
{
  // GaussianProcess class
  py::class_<GaussianProcess, std::shared_ptr<GaussianProcess>>
    gaussProc(m, "GaussianProcess");
  gaussProc
    .def(py::init<MeanFunctionBase&, KernelBase&>())
    .def(py::init<std::shared_ptr<MeanFunctionBase>, std::shared_ptr<KernelBase>>())
    .def("Condition", (GaussianProcess& (GaussianProcess::*)
                       (Eigen::Ref<const Eigen::MatrixXd> const&,
                        Eigen::Ref<const Eigen::MatrixXd> const&))
                        &GaussianProcess::Condition)
    .def("Condition", (GaussianProcess& (GaussianProcess::*)
                       (Eigen::Ref<const Eigen::MatrixXd> const&,
                        Eigen::Ref<const Eigen::MatrixXd> const&,
                        double)) &GaussianProcess::Condition)
    .def("Condition", (GaussianProcess& (GaussianProcess::*)
                       (std::shared_ptr<ObservationInformation>))
                       &GaussianProcess::Condition)
    .def("Optimize", &GaussianProcess::Optimize)
    .def("Predict", &GaussianProcess::Predict)
    .def("PredictMean", &GaussianProcess::PredictMean)
    .def("Sample", &GaussianProcess::Sample)
    .def("LogLikelihood", &GaussianProcess::LogLikelihood)
    .def("MarginalLogLikelihood", (double (GaussianProcess::*)())
     &GaussianProcess::MarginalLogLikelihood)
    .def("MarginalLogLikelihood", (double (GaussianProcess::*)
      (Eigen::Ref<Eigen::VectorXd>)) &GaussianProcess::MarginalLogLikelihood)
    .def("MarginalLogLikelihood", (double (GaussianProcess::*)
      (Eigen::Ref<Eigen::VectorXd>, bool)) &GaussianProcess::MarginalLogLikelihood)
    .def("Mean", &GaussianProcess::Mean)
    .def("Kernel", &GaussianProcess::Kernel)
    .def("Discretize", &GaussianProcess::Discretize);

  py::enum_<GaussianProcess::CovarianceType>(gaussProc, "CovarianceType")
        .value("DiagonalCov", GaussianProcess::CovarianceType::DiagonalCov)
        .value("BlockCov", GaussianProcess::CovarianceType::BlockCov)
        .value("FullCov", GaussianProcess::CovarianceType::FullCov)
        .value("NoCov", GaussianProcess::CovarianceType::NoCov)
        .export_values();

  // MeanFunctionBase class
  py::class_<MeanFunctionBase,
             std::shared_ptr<MeanFunctionBase>>
    meanFuncBase(m, "MeanFunctionBase");
  meanFuncBase
    //.def(py::init<unsigned, unsigned>())
    .def("Evaluate", &MeanFunctionBase::Evaluate)
    .def("Clone", &MeanFunctionBase::Clone)
    .def("GetPtr", &MeanFunctionBase::GetPtr)
    .def_readonly("inputDim", &MeanFunctionBase::inputDim)
    .def_readonly("coDim", &MeanFunctionBase::coDim);

  // ZeroMean class
  py::class_<ZeroMean, MeanFunctionBase, std::shared_ptr<ZeroMean>>
    zeroMean(m, "ZeroMean");
  zeroMean
    .def(py::init<unsigned, unsigned>())
    .def("Clone", &ZeroMean::Clone)
    .def("Evaluate", &ZeroMean::Evaluate);

  // LinearMean class
  py::class_<LinearMean, MeanFunctionBase, std::shared_ptr<LinearMean>>
    linMean(m, "LinearMean");
  linMean
    .def(py::init<double, double>())
    .def(py::init<Eigen::MatrixXd const&, Eigen::VectorXd const&>())
    .def("Clone", &LinearMean::Clone)
    .def("Evaluate", &LinearMean::Evaluate);

  /*// LinearTransformMean class
  py::class_<LinearTransformMean< >, MeanFunctionBase,
             std::shared_ptr<LinearTransformMean< >>>
    linTransMean(m, "LinearTransformMean");
  linTransMean
    .def(py::init<LinearOperator const&, MeanType const&>())
    .def("Clone", &LinearTransformMean< >::Clone)
    .def("Evaluate", &LinearTransformMean< >::Evaluate);*/

  /*// SumMean class
  py::class_<SumMean, MeanFunctionBase, std::shared_ptr<SumMean>>
    sumMean(m, "SumMean");
  sumMean
    .def(py::init<MeanType1 const&, MeanType2 const&>())
    .def("Clone", &SumMean::Clone)
    .def("Evaluate", &SumMean::Evaluate);*/

  // ObservationInformation class
  py::class_<ObservationInformation, std::shared_ptr<ObservationInformation>>
    obsInfo(m, "ObservationInformation");
  obsInfo
    .def(py::init<std::shared_ptr<muq::Modeling::LinearOperator>,
                           Eigen::Ref<const Eigen::VectorXd> const&,
                           Eigen::Ref<const Eigen::VectorXd> const&,
                           Eigen::Ref<const Eigen::MatrixXd> const&>());

    // virtual void FillSelfCov(std::shared_ptr<KernelBase> kernel,
    //                          Eigen::Ref<Eigen::MatrixXd> covBlock);
    //
    // virtual void FillCrossCov(Eigen::Ref<const Eigen::VectorXd> const& otherLoc,
    //                           std::shared_ptr<KernelBase>              kernel,
    //                           Eigen::Ref<Eigen::MatrixXd>              covBlock);
    //
    // virtual void FillCrossCov(std::shared_ptr<ObservationInformation> otherObs,
    //                           std::shared_ptr<KernelBase>             kernel,
    //                           Eigen::Ref<Eigen::MatrixXd>             covBlock);

  // DerivativeObservation class
  py::class_<DerivativeObservation, ObservationInformation,
             std::shared_ptr<DerivativeObservation>>
    derivObs(m, "DerivativeObservation");
  derivObs
    .def(py::init<std::shared_ptr<muq::Modeling::LinearOperator>,
                          Eigen::Ref<const Eigen::VectorXd> const&,
                          Eigen::Ref<const Eigen::VectorXd> const&,
                          Eigen::Ref<const Eigen::MatrixXd> const&,
                          std::vector<std::vector<int>>>());

  //StateSpaceGP class
  py::class_<StateSpaceGP, GaussianProcess, std::shared_ptr<StateSpaceGP>>
    stateSpace(m, "StateSpaceGP");
  stateSpace
    .def(py::init<MeanFunctionBase&, KernelBase&>())
    .def(py::init<MeanFunctionBase&, KernelBase&, boost::property_tree::ptree>())
    .def(py::init<std::shared_ptr<MeanFunctionBase>, std::shared_ptr<KernelBase>>())
    .def(py::init<std::shared_ptr<MeanFunctionBase>, std::shared_ptr<KernelBase>,
                  boost::property_tree::ptree>())
    .def("Sample", &StateSpaceGP::Sample)
    .def("Predict", &StateSpaceGP::Predict)
    .def("PredictMean", &StateSpaceGP::PredictMean)
    .def("LogLikelihood", &StateSpaceGP::LogLikelihood)
    .def("MarginalLogLikelihood", &StateSpaceGP::MarginalLogLikelihood)
    .def("GetSDE", &StateSpaceGP::GetSDE)
    .def("GetObs", &StateSpaceGP::GetObs)
    .def("SetObs", &StateSpaceGP::SetObs)
    .def("GetCov", &StateSpaceGP::GetCov)
    .def_readonly("stateDim", &StateSpaceGP::stateDim);
    //.def("Concatenate", &StateSpaceGP::Concatenate);
}
