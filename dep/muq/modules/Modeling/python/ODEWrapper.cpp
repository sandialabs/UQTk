#include "AllClassWrappers.h"

#include "MUQ/Utilities/PyDictConversion.h"

#include "MUQ/Modeling/ODE.h"

using namespace muq::Modeling;
using namespace muq::Utilities;
namespace py = pybind11;

void PythonBindings::ODEWrapper(py::module &m) {
  py::class_<ODE, ModPiece, std::shared_ptr<ODE> > ode(m, "ODE");
  ode.def(py::init( [] (std::shared_ptr<ModPiece> const& rhs, py::dict const& d) { return new ODE(rhs, ConvertDictToPtree(d)); }));
#if MUQ_HAS_PARCER==1
  ode.def(py::init( [] (std::shared_ptr<ModPiece> const& rhs, py::dict const& d, std::shared_ptr<parcer::Communicator> const& comm) { return new ODE(rhs, ConvertDictToPtree(d), comm); }));
#endif
}
