#ifndef UTILITIES_ALLCLASSWRAPPERS_H_
#define UTILITIES_ALLCLASSWRAPPERS_H_

#include <pybind11/pybind11.h>

namespace muq{
  namespace Utilities{

    namespace PythonBindings{

      //void HDF5Wrapper(pybind11::module &m);
      //void LinearAlgebraWrapper(pybind11::module &m);
      void MultiIndicesWrapper(pybind11::module &m);
      void GeneralUtilitiesWrapper(pybind11::module &m);
      
    }
  }
}


#endif