#ifndef MODELING_ALLCLASSWRAPPERS_H_
#define MODELING_ALLCLASSWRAPPERS_H_

#include <pybind11/pybind11.h>

namespace muq{
  namespace Modeling{
    namespace PythonBindings{

      void WorkPieceWrapper(pybind11::module &m);
      void ModPieceWrapper(pybind11::module &m);
      void DistributionWrapper(pybind11::module &m);
      void CwiseUnaryOperatorsWrapper(pybind11::module &m);
      void LinearOperatorWrapper(pybind11::module &m);
      void ODEWrapper(pybind11::module &m);
    }
  }
}


#endif
