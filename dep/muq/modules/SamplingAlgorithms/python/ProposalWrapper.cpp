#include "AllClassWrappers.h"

#include "MUQ/SamplingAlgorithms/AMProposal.h"
#include "MUQ/SamplingAlgorithms/MCMCProposal.h"
#include "MUQ/SamplingAlgorithms/MHProposal.h"
#include "MUQ/SamplingAlgorithms/CrankNicolsonProposal.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

//using namespace muq::SamplingAlgorithms::PythonBindings;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;
namespace py = pybind11;

void PythonBindings::ProposalWrapper(py::module &m) {
  py::class_<MCMCProposal, std::shared_ptr<MCMCProposal>> mcmcPro(m, "MCMCProposal");
  mcmcPro
    .def("Sample", &MCMCProposal::Sample)
    .def("LogDensity", &MCMCProposal::LogDensity);
    //.def("Construct", &MCMCProposal::Construct)
    //.def("GetMCMCProposalMap", &MCMCProposal::GetMCMCProposalMap);

  py::class_<MHProposal, MCMCProposal, std::shared_ptr<MHProposal>> mhPro(m, "MHProposal");
  mhPro
    .def(py::init([](py::dict d, std::shared_ptr<AbstractSamplingProblem> prob) {return new MHProposal(ConvertDictToPtree(d), prob);} ))
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> prob, std::shared_ptr<muq::Modeling::GaussianBase> gauss) { return new MHProposal(ConvertDictToPtree(d), prob, gauss);}));


  py::class_<CrankNicolsonProposal, MCMCProposal, std::shared_ptr<CrankNicolsonProposal>> cnPro(m, "CrankNicolsonProposal");
  cnPro
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> prob) {return new CrankNicolsonProposal(ConvertDictToPtree(d), prob);} ))
    .def(py::init( [](py::dict d, std::shared_ptr<AbstractSamplingProblem> prob, std::shared_ptr<muq::Modeling::GaussianBase> gauss) { return new CrankNicolsonProposal(ConvertDictToPtree(d), prob, gauss);}));


  py::class_<AMProposal, MHProposal, std::shared_ptr<AMProposal>> amPro(m, "AMProposal");
  amPro
    .def(py::init<boost::property_tree::ptree const&, std::shared_ptr<AbstractSamplingProblem>>())
    .def("Adapt", &AMProposal::Adapt);

}
