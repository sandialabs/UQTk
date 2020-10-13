#include "AllClassWrappers.h"

#include "MUQ/Utilities/HDF5/Attributes.h"
#include "MUQ/Utilities/HDF5/BlockDataset.h"
#include "MUQ/Utilities/HDF5/H5Object.h"
#include "MUQ/Utilities/HDF5/HDF5File.h"
#include "MUQ/Utilities/HDF5/HDF5Types.h"
#include "MUQ/Utilities/HDF5/Pathtools.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::HDF5Wrapper(py::module &m)
{
    
}