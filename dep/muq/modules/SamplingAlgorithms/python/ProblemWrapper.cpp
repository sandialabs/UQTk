#include "AllClassWrappers.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingProblem.h"
#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"
#include "MUQ/SamplingAlgorithms/SamplingState.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Utilities;
using namespace muq::SamplingAlgorithms;
namespace py = pybind11;

void PythonBindings::ProblemWrapper(py::module &m) {
  py::class_<AbstractSamplingProblem, std::shared_ptr<AbstractSamplingProblem>> absSamp(m, "AbstractSamplingProblem");
  absSamp
    .def("LogDensity", &AbstractSamplingProblem::LogDensity)
    .def_readonly("numBlocks", &AbstractSamplingProblem::numBlocks)
    .def_readonly("blockSizes", &AbstractSamplingProblem::blockSizes);

  py::class_<SamplingProblem, AbstractSamplingProblem, std::shared_ptr<SamplingProblem>> sampProb(m, "SamplingProblem");
  sampProb
    .def(py::init<std::shared_ptr<muq::Modeling::ModPiece>>())
    .def("LogDensity", &SamplingProblem::LogDensity)
    .def("GradLogDensity", &SamplingProblem::GradLogDensity)
    .def("GetDistribution", &SamplingProblem::GetDistribution);

    py::class_<ExpensiveSamplingProblem, SamplingProblem, std::shared_ptr<ExpensiveSamplingProblem>> expenProb(m, "ExpensiveSamplingProblem");
    expenProb
      .def(py::init( [] (std::shared_ptr<muq::Modeling::ModPiece> target, py::dict d) { return new ExpensiveSamplingProblem(target, ConvertDictToPtree(d)); }))
      .def(py::init( [] (std::shared_ptr<muq::Modeling::ModPiece> target, Eigen::VectorXd const& centroid, py::dict d) { return new ExpensiveSamplingProblem(target, centroid, ConvertDictToPtree(d)); }))
      .def("LogDensity", &ExpensiveSamplingProblem::LogDensity)
      .def("CacheSize", &ExpensiveSamplingProblem::CacheSize);

#if MUQ_HAS_PARCER==1
      expenProb.def(py::init( [] (std::shared_ptr<muq::Modeling::ModPiece> target, Eigen::VectorXd const& centroid, py::dict d, std::shared_ptr<parcer::Communicator> comm) { return new ExpensiveSamplingProblem(target, centroid, ConvertDictToPtree(d), comm); }));
#endif
}
