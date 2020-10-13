#ifndef SAMPLINGALGORITHMS_ALLCLASSWRAPPERS_H_
#define SAMPLINGALGORITHMS_ALLCLASSWRAPPERS_H_

#include <pybind11/pybind11.h>

namespace muq{
  namespace SamplingAlgorithms{
    namespace PythonBindings{

      void KernelWrapper(pybind11::module &m);
      void ProposalWrapper(pybind11::module &m);
      void SampleWrapper(pybind11::module &m);
      void MCMCWrapper(pybind11::module &m);
      void ProblemWrapper(pybind11::module &m);
    }
  }
}


#endif
