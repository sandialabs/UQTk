#include "AllClassWrappers.h"

#include "MUQ/Utilities/RandomGenerator.h"
#include "MUQ/Utilities/StringUtilities.h"
#include "MUQ/Utilities/WaitBar.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Utilities::PythonBindings;
namespace py = pybind11;

void muq::Utilities::PythonBindings::GeneralUtilitiesWrapper(py::module &m)
{
  py::class_<RandomGenerator, std::shared_ptr<RandomGenerator>> randGen(m, "RandomGenerator");
  randGen
  //.def_static("GetNormal", (double (RandomGenerator::*)()) &RandomGenerator::GetNormal)
  .def_static("GetNormal", (double (*)()) &RandomGenerator::GetNormal)
  .def_static("GetNormal", (Eigen::MatrixXd (*)(int)) &RandomGenerator::GetNormal)
  .def_static("GetNormal", (Eigen::MatrixXd (*)(int, int)) &RandomGenerator::GetNormal)
  
  .def_static("GetUniform", (double (*)()) &RandomGenerator::GetUniform)
  .def_static("GetUniform", (Eigen::MatrixXd (*)(int)) &RandomGenerator::GetUniform)
  .def_static("GetUniform", (Eigen::MatrixXd (*)(int, int)) &RandomGenerator::GetUniform)
  
  .def_static("GetGamma", (double (*)(double const, double const)) &RandomGenerator::GetGamma)
  .def_static("GetGamma", (Eigen::MatrixXd (*)(double const, double const, int)) &RandomGenerator::GetGamma)
  .def_static("GetGamma", (Eigen::MatrixXd (*)(double const, double const, int, int)) &RandomGenerator::GetGamma)
  
  .def_static("GetUniformInt", (int (*)(int, int)) &RandomGenerator::GetUniformInt)
  .def_static("GetUniformInt", (Eigen::MatrixXi (*)(int, int, int)) &RandomGenerator::GetUniformInt)
  .def_static("GetUniformInt", (Eigen::MatrixXi (*)(int, int, int, bool)) &RandomGenerator::GetUniformInt)
  .def_static("GetUniformInt", (Eigen::MatrixXi (*)(int, int, int, int)) &RandomGenerator::GetUniformInt)
  .def_static("GetUniformInt", (Eigen::MatrixXi (*)(int, int, int, int, bool)) &RandomGenerator::GetUniformInt)
  
  .def_static("SetSeed", &RandomGenerator::SetSeed)
  .def_static("CopyGenerator", &RandomGenerator::CopyGenerator)
  .def_static("SetGenerator", &RandomGenerator::SetGenerator);
}