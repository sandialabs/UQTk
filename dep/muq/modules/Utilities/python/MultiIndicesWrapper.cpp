#include "AllClassWrappers.h"

#include "MUQ/Utilities/MultiIndices/MultiIndex.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexFactory.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexLimiter.h"
#include "MUQ/Utilities/MultiIndices/MultiIndexSet.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <string>

#include <functional>
#include <vector>

using namespace muq::Utilities::PythonBindings;
namespace py = pybind11;

void muq::Utilities::PythonBindings::MultiIndicesWrapper(py::module &m)
{
  py::class_<MultiIndex, std::shared_ptr<MultiIndex>> multiI(m, "MultiIndex");
  multiI
    .def(py::init<unsigned>())
    .def(py::init<Eigen::RowVectorXi const&>())
    .def(py::init<std::initializer_list<unsigned> const&>())
    .def_static("Copy", &MultiIndex::Copy)
    .def("GetVector", &MultiIndex::GetVector)
    .def("Sum", &MultiIndex::Sum)
    .def("Max", &MultiIndex::Max)
    .def("SetValue", &MultiIndex::SetValue)
    .def("GetValue", &MultiIndex::GetValue)
    .def("SetLength", &MultiIndex::SetLength)
    .def("GetLength", &MultiIndex::GetLength)
    .def("GetNzBegin", &MultiIndex::GetNzBegin)
    .def("GetNzEnd", &MultiIndex::GetNzEnd);

  py::class_<MultiIndexLimiter, std::shared_ptr<MultiIndexLimiter>> multiILim(m, "MultiIndexLimiter");
  multiILim
    .def("IsFeasible", &MultiIndexLimiter::IsFeasible);

  py::class_<TotalOrderLimiter, MultiIndexLimiter, std::shared_ptr<TotalOrderLimiter>> totalOrd(m, "TotalOrderLimiter");
  totalOrd
    .def(py::init<unsigned int>());

  py::class_<DimensionLimiter, MultiIndexLimiter, std::shared_ptr<DimensionLimiter>> dimLim(m, "DimensionLimiter");
  dimLim
    .def(py::init<unsigned int, unsigned int>())
    .def("IsFeasible", &DimensionLimiter::IsFeasible);

  py::class_<MaxOrderLimiter, MultiIndexLimiter, std::shared_ptr<MaxOrderLimiter>> maxLim(m, "MaxOrderLimiter");
  maxLim
    .def(py::init<unsigned int>())
    .def(py::init<Eigen::VectorXi const&>())
    .def("IsFeasible", &MaxOrderLimiter::IsFeasible);

  py::class_<NoLimiter, MultiIndexLimiter, std::shared_ptr<NoLimiter>> noLim(m, "NoLimiter");
  noLim
    .def(py::init<>())
    .def("IsFeasible", &NoLimiter::IsFeasible);

  py::class_<AndLimiter, MultiIndexLimiter, std::shared_ptr<AndLimiter>> andLim(m, "AndLimiter");
  andLim
    .def(py::init<std::shared_ptr<MultiIndexLimiter>, std::shared_ptr<MultiIndexLimiter>>())
    .def("IsFeasible", &AndLimiter::IsFeasible);

  py::class_<OrLimiter, MultiIndexLimiter, std::shared_ptr<OrLimiter>> orLim(m, "OrLimiter");
  orLim
    .def(py::init<std::shared_ptr<MultiIndexLimiter>, std::shared_ptr<MultiIndexLimiter>>())
    .def("IsFeasible", &OrLimiter::IsFeasible);

  py::class_<XorLimiter, MultiIndexLimiter, std::shared_ptr<XorLimiter>> xorLim(m, "XorLimiter");
  xorLim
    .def(py::init<std::shared_ptr<MultiIndexLimiter>, std::shared_ptr<MultiIndexLimiter>>())
    .def("IsFeasible", &XorLimiter::IsFeasible);

  py::class_<MultiIndexFactory, std::shared_ptr<MultiIndexFactory>> multiIFac(m, "MultiIndexFactory");
  multiIFac
    .def_static("CreateTotalOrder", &MultiIndexFactory::CreateTotalOrder, py::arg("length"), py::arg("maxOrder"), py::arg("minOrder")=0, py::arg("limiter")=std::make_shared<NoLimiter>())
    .def_static("CreateTriTotalOrder", &MultiIndexFactory::CreateTriTotalOrder, py::arg("length"), py::arg("maxOrder"), py::arg("minOrder")=0, py::arg("limiter")=std::make_shared<NoLimiter>())
    .def_static("CreateHyperbolic", &MultiIndexFactory::CreateHyperbolic, py::arg("length"), py::arg("maxOrder"), py::arg("q")=0, py::arg("limiter")=std::make_shared<NoLimiter>())
    .def_static("CreateTriHyperbolic", &MultiIndexFactory::CreateTriHyperbolic, py::arg("length"), py::arg("maxOrder"), py::arg("q")=0, py::arg("limiter")=std::make_shared<NoLimiter>())
    //.def("CentralMoment", &SampleCollection::CentralMoment, py::arg("order"), py::arg("blockDim") = -1)
    //.def_static("CreateFullTensor", (std::shared_ptr<MultiIndexSet> (MultiIndexFactory::*)(unsigned int const, unsigned int const,  std::shared_ptr<MultiIndexLimiter>)) &MultiIndexFactory::CreateFullTensor)
    //.def_static("CreateFullTensor", (std::shared_ptr<MultiIndexSet> (MultiIndexFactory::*)(const Eigen::RowVectorXi&, std::shared_ptr<MultiIndexLimiter>)) &MultiIndexFactory::CreateFullTensor)
    .def_static("CreateTriHyperbolic", &MultiIndexFactory::CreateTriHyperbolic);

  py::class_<MultiIndexSet,std::shared_ptr<MultiIndexSet>> multiSet(m, "MultiIndexSet");
  multiSet
    .def(py::init<const unsigned , std::shared_ptr<MultiIndexLimiter>>())
    .def_static("CloneExisting", &MultiIndexSet::CloneExisting)
    .def("SetLimiter", &MultiIndexSet::SetLimiter)
    .def("GetLimiter", &MultiIndexSet::GetLimiter)
    .def("IndexToMulti", &MultiIndexSet::IndexToMulti)
    .def("MultiToIndex", &MultiIndexSet::MultiToIndex)
    .def("GetMultiLength", &MultiIndexSet::GetMultiLength)
    .def("GetMaxOrders", &MultiIndexSet::GetMaxOrders)
    .def("at", &MultiIndexSet::at)
    .def("Size", &MultiIndexSet::Size)
    .def("Union", &MultiIndexSet::Union)
    .def("Activate", (void (MultiIndexSet::*)(std::shared_ptr<MultiIndex> const&)) &MultiIndexSet::Activate)
    .def("AddActive", &MultiIndexSet::AddActive)
    .def("Expand", &MultiIndexSet::Expand)
    .def("ForciblyExpand", &MultiIndexSet::ForciblyExpand)
    .def("ForciblyActivate", (std::vector<unsigned> (MultiIndexSet::*)(std::shared_ptr<MultiIndex> const&)) &MultiIndexSet::ForciblyActivate)
    .def("GetAdmissibleForwardNeighbors", &MultiIndexSet::GetAdmissibleForwardNeighbors)
    .def("IsAdmissible", (bool (MultiIndexSet::*)(std::shared_ptr<MultiIndex> const&) const) &MultiIndexSet::IsAdmissible)
    .def("IsExpandable", &MultiIndexSet::IsExpandable)
    .def("IsActive", (bool (MultiIndexSet::*)(std::shared_ptr<MultiIndex> const&) const) &MultiIndexSet::IsActive);

}
