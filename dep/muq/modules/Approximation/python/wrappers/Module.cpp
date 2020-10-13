#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "AllClassWrappers.h"

using namespace muq::Approximation::PythonBindings;
namespace py = pybind11;

PYBIND11_MODULE(pymuqApproximationWrappers, m) {

    KernelWrapper(m);
    GaussianWrapper(m);
    PolynomialsWrapper(m);
    KLWrapper(m);
    QuadratureWrapper(m);
    PolynomialChaosWrapper(m);
}
