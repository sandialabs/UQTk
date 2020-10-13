#include "AllClassWrappers.h"

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"
#include "MUQ/Approximation/GaussianProcesses/KernelImpl.h"

#include "MUQ/Approximation/GaussianProcesses/ConcatenateKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ConstantKernel.h"
#include "MUQ/Approximation/GaussianProcesses/LinearTransformKernel.h"
#include "MUQ/Approximation/GaussianProcesses/MaternKernel.h"
#include "MUQ/Approximation/GaussianProcesses/PeriodicKernel.h"
#include "MUQ/Approximation/GaussianProcesses/ProductKernel.h"
#include "MUQ/Approximation/GaussianProcesses/SquaredExpKernel.h"
#include "MUQ/Approximation/GaussianProcesses/SumKernel.h"
#include "MUQ/Approximation/GaussianProcesses/WhiteNoiseKernel.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Approximation::PythonBindings;
using namespace muq::Approximation;
namespace py = pybind11;

void muq::Approximation::PythonBindings::KernelWrapper(py::module &m)
{
    // KernelBase class
    py::class_<KernelBase,
               std::shared_ptr<KernelBase>> kernBase(m, "CovarianceKernelBase");
    kernBase
      //.def(py::init<unsigned, unsigned, unsigned>())
      //.def(py::init<unsigned, std::vector<unsigned>, unsigned, unsigned>())
      .def("GetSeperableComponents", &KernelBase::GetSeperableComponents)
      .def("Evaluate", &KernelBase::Evaluate)
      .def("BuildCovariance", (Eigen::MatrixXd (KernelBase::*)
                              (Eigen::MatrixXd const&) const)
                              &KernelBase::BuildCovariance)
      .def("BuildCovariance", (Eigen::MatrixXd (KernelBase::*)
                              (Eigen::MatrixXd const&,
                               Eigen::MatrixXd const&) const)
                              &KernelBase::BuildCovariance)
      .def("FillCovariance", (void (KernelBase::*)
                             (Eigen::MatrixXd const&,
                              Eigen::MatrixXd const&,
                              Eigen::Ref<Eigen::MatrixXd>) const)
                             &KernelBase::FillCovariance)
      .def("FillCovariance", (void (KernelBase::*)
                             (Eigen::MatrixXd const&,
                              Eigen::Ref<Eigen::MatrixXd>) const)
                             &KernelBase::FillCovariance)
      .def("FillDerivCovariance", &KernelBase::FillDerivCovariance)
      .def("GetPosDerivative", &KernelBase::GetPosDerivative)
      .def("GetParamBounds", &KernelBase::GetParamBounds)
      .def("GetParams", &KernelBase::GetParams)
      .def("SetParams", &KernelBase::SetParams)
      .def("Clone", &KernelBase::Clone)
      .def("GetStateSpace", &KernelBase::GetStateSpace)
      .def_readonly("dimInds", &KernelBase::dimInds)
      .def_readonly("inputDim", &KernelBase::inputDim)
      .def_readonly("coDim", &KernelBase::coDim)
      .def_readonly("numParams", &KernelBase::numParams)
      .def("FillPosDerivBlock", &KernelBase::FillPosDerivBlock)
      .def("FillBlock", &KernelBase::FillBlock)
      .def("__mul__", [](std::shared_ptr<KernelBase> a, std::shared_ptr<KernelBase> b) {return a * b;}, py::is_operator())
      .def("__add__", [](std::shared_ptr<KernelBase> a, std::shared_ptr<KernelBase> b) {return a + b;}, py::is_operator());

    // ConcatenateKernel class
    py::class_<ConcatenateKernel, KernelBase, std::shared_ptr<ConcatenateKernel>>
      concatKern(m, "ConcatenateKernel");
    concatKern
      .def(py::init<std::vector<std::shared_ptr<KernelBase>> const&>())
      .def(py::init<std::shared_ptr<KernelBase> const&,
                    std::shared_ptr<KernelBase> const&>())
      .def("Clone", &ConcatenateKernel::Clone)
      .def("FillBlock", &ConcatenateKernel::FillBlock)
      .def("FillPosDerivBlock", &ConcatenateKernel::FillPosDerivBlock);

    // KernelImpl<ConstantKernel> class
    py::class_<KernelImpl<ConstantKernel>, KernelBase,
               std::shared_ptr<KernelImpl<ConstantKernel>>>
      kernImplConst(m, "ConstantKernelImpl");
    kernImplConst
      .def(py::init<unsigned, unsigned, unsigned>())
      .def(py::init<unsigned, std::vector<unsigned>, unsigned, unsigned>())
      .def("Clone", &KernelImpl<ConstantKernel>::Clone)
      .def("FillBlock", &KernelImpl<ConstantKernel>::FillBlock)
      .def("FillPosDerivBlock", &KernelImpl<ConstantKernel>::FillPosDerivBlock)
      .def("FillPosDerivBlockImpl", &KernelImpl<ConstantKernel>::FillPosDerivBlockImpl);

    // ConstantKernel class
    py::class_<ConstantKernel, KernelImpl<ConstantKernel>, std::shared_ptr<ConstantKernel>>
      constKern(m, "ConstantKernel");
    constKern
      .def(py::init<unsigned, const double>())
      .def(py::init<unsigned, const double, const Eigen::Vector2d>())
      .def(py::init<unsigned, std::vector<unsigned>, const double>())
      .def(py::init<unsigned, std::vector<unsigned>, const double,
                    const Eigen::Vector2d>())
      .def(py::init<unsigned, Eigen::MatrixXd const&>())
      .def(py::init<unsigned, Eigen::MatrixXd const&, const Eigen::Vector2d>())
      .def(py::init<unsigned, std::vector<unsigned>, Eigen::MatrixXd const&>())
      .def(py::init<unsigned, std::vector<unsigned>, Eigen::MatrixXd const&,
                    const Eigen::Vector2d>());


    // LinearTransformKernel class
    py::class_<LinearTransformKernel, KernelBase,
               std::shared_ptr<LinearTransformKernel>>
      linTransKern(m, "LinearTransformKernel");
    linTransKern
      .def(py::init<Eigen::MatrixXd const&, std::shared_ptr<KernelBase>>())
      .def("FillBlock", &LinearTransformKernel::FillBlock)
      .def("FillPosDerivBlock", &LinearTransformKernel::FillPosDerivBlock)
      .def("Clone", &LinearTransformKernel::Clone);

    // KernelImpl<MaternKernel> class
    py::class_<KernelImpl<MaternKernel>, KernelBase,
               std::shared_ptr<KernelImpl<MaternKernel>>>
      kernImplMat(m, "MaternKernelImpl");
    kernImplMat
      .def(py::init<unsigned, unsigned, unsigned>())
      .def(py::init<unsigned, std::vector<unsigned>, unsigned, unsigned>())
      .def("Clone", &KernelImpl<MaternKernel>::Clone)
      .def("FillBlock", &KernelImpl<MaternKernel>::FillBlock)
      .def("FillPosDerivBlock", &KernelImpl<MaternKernel>::FillPosDerivBlock)
      .def("FillPosDerivBlockImpl", &KernelImpl<MaternKernel>::FillPosDerivBlockImpl);

    // MaternKernel class
    py::class_<MaternKernel, KernelImpl<MaternKernel>,
               std::shared_ptr<MaternKernel>>
      matKern(m, "MaternKernel");
    matKern
      .def(py::init<unsigned, std::vector<unsigned>, double, double, double>())
      .def(py::init<unsigned, std::vector<unsigned>, double, double, double,Eigen::Vector2d, Eigen::Vector2d>())
      .def(py::init<unsigned, double, double, double>())
      .def(py::init<unsigned, double, double, double, Eigen::Vector2d,Eigen::Vector2d>())
      .def("GetStateSpace", &MaternKernel::GetStateSpace);

    // KernelImpl<PeriodicKernel> class
    py::class_<KernelImpl<PeriodicKernel>, KernelBase,
               std::shared_ptr<KernelImpl<PeriodicKernel>>>
      kernImplPer(m, "PeriodicKernelImpl");
    kernImplPer
      .def(py::init<unsigned, unsigned, unsigned>())
      .def(py::init<unsigned, std::vector<unsigned>, unsigned, unsigned>())
      .def("Clone", &KernelImpl<PeriodicKernel>::Clone)
      .def("FillBlock", &KernelImpl<PeriodicKernel>::FillBlock)
      .def("FillPosDerivBlock", &KernelImpl<PeriodicKernel>::FillPosDerivBlock)
      .def("FillPosDerivBlockImpl", &KernelImpl<PeriodicKernel>::FillPosDerivBlockImpl);

    // PeriodicKernel class
    py::class_<PeriodicKernel, KernelImpl<PeriodicKernel>,
               std::shared_ptr<PeriodicKernel>>
      periodKern(m, "PeriodicKernel");
    periodKern
      .def(py::init<unsigned, std::vector<unsigned>, double, double, double>())
      .def(py::init<unsigned, std::vector<unsigned>, double, double, double, Eigen::Vector2d, Eigen::Vector2d, Eigen::Vector2d>())
      .def(py::init<unsigned, double, double, double>())
      .def(py::init<unsigned, double, double, double,Eigen::Vector2d, Eigen::Vector2d, Eigen::Vector2d>())
      //.def("FillBlockImpl", &PeriodicKernel::FillBlockImpl)
      .def("GetStateSpace", &PeriodicKernel::GetStateSpace);

    // ProductKernel class
    py::class_<ProductKernel, KernelBase, std::shared_ptr<ProductKernel>>
      productKern(m, "ProductKernel");
    productKern
      .def(py::init<std::shared_ptr<KernelBase>, std::shared_ptr<KernelBase>>())
      .def("FillBlock", &ProductKernel::FillBlock)
      .def("FillPosDerivBlock", &ProductKernel::FillPosDerivBlock)
      .def("Clone", &ProductKernel::Clone)
      .def("GetSeperableComponents", &ProductKernel::GetSeperableComponents)
      .def("GetStateSpace", &ProductKernel::GetStateSpace);

    // KernelImpl<SquaredExpKernel> class
    py::class_<KernelImpl<SquaredExpKernel>, KernelBase,
               std::shared_ptr<KernelImpl<SquaredExpKernel>>>
      kernImplSE(m, "SquaredExpKernelImpl");
    kernImplSE
      .def(py::init<unsigned, unsigned, unsigned>())
      .def(py::init<unsigned, std::vector<unsigned>, unsigned, unsigned>())
      .def("Clone", &KernelImpl<SquaredExpKernel>::Clone)
      .def("FillBlock", &KernelImpl<SquaredExpKernel>::FillBlock)
      .def("FillPosDerivBlock", &KernelImpl<SquaredExpKernel>::FillPosDerivBlock)
      .def("FillPosDerivBlockImpl", &KernelImpl<SquaredExpKernel>::FillPosDerivBlockImpl);

    // SquaredExpKernel class
    py::class_<SquaredExpKernel, KernelImpl<SquaredExpKernel>,
               std::shared_ptr<SquaredExpKernel>>
      sqExpKern(m, "SquaredExpKernel");
    sqExpKern
      .def(py::init<unsigned, std::vector<unsigned>, double, double>())
      .def(py::init<unsigned, std::vector<unsigned>, double, double, Eigen::Vector2d, Eigen::Vector2d>())
      .def(py::init<unsigned, double, double, Eigen::Vector2d, Eigen::Vector2d>())
      .def(py::init<unsigned, double, double>());
      //.def("FillBlockImpl", &SquaredExpKernel::FillBlockImpl);

    // SumKernel class
    py::class_<SumKernel, KernelBase, std::shared_ptr<SumKernel>>
      sumKern(m, "SumKernel");
    sumKern
      .def(py::init<std::shared_ptr<KernelBase>, std::shared_ptr<KernelBase>>())
      .def("Clone", &SumKernel::Clone)
      .def("FillBlock", &SumKernel::FillBlock)
      .def("FillPosDerivBlock", &SumKernel::FillPosDerivBlock)
      .def("GetStateSpace", &SumKernel::GetStateSpace);

    // KernelImpl<WhiteNoiseKernel> class
    py::class_<KernelImpl<WhiteNoiseKernel>, KernelBase,
               std::shared_ptr<KernelImpl<WhiteNoiseKernel>>>
      kernImplWN(m, "WhiteNoiseKernelImpl");
    kernImplWN
      .def(py::init<unsigned, unsigned, unsigned>())
      .def(py::init<unsigned, std::vector<unsigned>, unsigned, unsigned>())
      .def("Clone", &KernelImpl<WhiteNoiseKernel>::Clone)
      .def("FillBlock", &KernelImpl<WhiteNoiseKernel>::FillBlock)
      .def("FillPosDerivBlock", &KernelImpl<WhiteNoiseKernel>::FillPosDerivBlock)
      .def("FillPosDerivBlockImpl", &KernelImpl<WhiteNoiseKernel>::FillPosDerivBlockImpl);

    // WhiteNoiseKernel class
    py::class_<WhiteNoiseKernel, KernelImpl<WhiteNoiseKernel>,
               std::shared_ptr<WhiteNoiseKernel>>
      whiteNoiseKern(m, "WhiteNoiseKernel");
    whiteNoiseKern
      .def(py::init<unsigned, const double>())
      .def(py::init<unsigned, const double, const Eigen::Vector2d >());
      //.def("FillBlockImpl", &WhiteNoiseKernel::FillBlockImpl);
}
