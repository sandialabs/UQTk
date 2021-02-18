#include <pybind11/pybind11.h>
#include <>

PYBIND11_MODULE(uqtkarray, m) {
    m.class_<Array1D>(m,"Array1D")
        .def(py::init<int &,std::vector &>())
        .def("Clear",&Array1D::Clear)
        .def("XSize",&Array1D::XSize)
        .def("Length",&Array1D::Length)
        .def("Resize",static_cast<void (Array1D::*)(const int &)>(&Array1D::Resize))
        .def("Resize",static_cast<void (Array1D::*)(const int &,const T&)>(&Array1D::Resize))
        .def("SetValue",&Array1D::SetValue)
        .def("PushBack",&Array1D::PushBack)
        .def("GetArrayPointer",&Array1D::GetArrayPointer)
        .def("GetConstArrayPointer",&Array1D::GetConstArrayPointer)
        .def("insert",static_cast<void (Array1D::*)(Array1D<T>&,int)>(&Array1D::insert))
        .def("insert",static_cast<void (Array1D::*)(const T&,int)>(&Array1D::insert))
        .def("erase",&Array1D::erase)
        .def("DumpBinary",static_cast<void (Array1D::*)(FILE*)>(&Array1D::DumpBinary))
        .def("DumpBinary",static_cast<void (Array1D::*)(char*)>(&Array1D::DumpBinary))
        .def("ReadBinary",&Array1D::ReadBinary)
        .def(
}
