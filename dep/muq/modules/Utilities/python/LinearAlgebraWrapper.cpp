#include "AllClassWrappers.h"

#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"
#include "MUQ/Modeling/LinearAlgebra/BlockDiagonalOperator.h"
#include "MUQ/Modeling/LinearAlgebra/BlockRowOperator.h"
#include "MUQ/Modeling/LinearAlgebra/CompanionMatrix.h"
#include "MUQ/Modeling/LinearAlgebra/ConcatenateOperator.h"
#include "MUQ/Modeling/LinearAlgebra/DiagonalOperator.h"
#include "MUQ/Modeling/LinearAlgebra/EigenMatrixAlgebra.h"
#include "MUQ/Modeling/LinearAlgebra/EigenVectorAlgebra.h"
#include "MUQ/Modeling/LinearAlgebra/KroneckerProductOperator.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"
#include "MUQ/Modeling/LinearAlgebra/ProductOperator.h"
#include "MUQ/Modeling/LinearAlgebra/ScalarAlgebra.h"
#include "MUQ/Modeling/LinearAlgebra/SumOperator.h"
#include "MUQ/Modeling/LinearAlgebra/SundialsAlgebra.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::SamplingAlgorithms::PythonBindings;
namespace py = pybind11;

void muq::SamplingAlgorithms::PythonBindings::LinearAlgebraWrapper(py::module &m)
{
    
}