#ifndef APPROXIMATION_ALLCLASSWRAPPERS_H_
#define APPROXIMATION_ALLCLASSWRAPPERS_H_

#include <pybind11/pybind11.h>

namespace muq{
  namespace Approximation{
    namespace PythonBindings{

      void KernelWrapper(pybind11::module &m);
      void GaussianWrapper(pybind11::module &m);
      void PolynomialsWrapper(pybind11::module &m);
      void KLWrapper(pybind11::module &m);
      void QuadratureWrapper(pybind11::module &m);
      void PolynomialChaosWrapper(pybind11::module &m);

    }
  }
}


#endif
