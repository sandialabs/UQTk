#ifndef OPTIMIZATION_ALLCLASSWRAPPERS_H_
#define OPTIMIZATION_ALLCLASSWRAPPERS_H_

#include <pybind11/pybind11.h>

namespace muq {
  namespace Optimization {
    namespace PythonBindings {
      void CostFunctionWrapper(pybind11::module &m);
      void OptimizationWrapper(pybind11::module &m);
    } // namespace PythonBindings
  } // namespace Optimization
} // namespace muq

#endif
