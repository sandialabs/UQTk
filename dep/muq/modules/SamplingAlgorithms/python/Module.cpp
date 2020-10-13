#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

PYBIND11_MODULE(pymuqSamplingAlgorithms, m) {

    KernelWrapper(m);
    ProposalWrapper(m);
    SampleWrapper(m);
    MCMCWrapper(m);
    ProblemWrapper(m);
}
