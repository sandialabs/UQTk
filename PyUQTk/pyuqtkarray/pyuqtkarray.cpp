#include <pybind11/pybind11.h>

#include "Array1D.h"
#include "Array2D.h"
#include "arrayio.h"
#include "arraytool.h"

using Array1D<int> = Array1D::Array1D<int>;
using Array1D<double> = Array1D::Array1D<double>;
using Array2D<int> = Array2D::Array2D<int>;
using Array2D<double> = Array2D::Array2D<double>;

namespace py=pybind11;

PYBIND11_MODULE(uqtkarray, m) {
    py::class_<Array1D<int>>(m, "Array1D<int>")
      .def(py::init<>())
      .def(py::init<const int&>())
      .def(py::init<const int&,const int&>())
      .def("Assign", &Array1D::operator=)
      .def(py::init<const Array1D &>())
      .def("Clear",&Array1D::Clear)
      .def("XSize",&Array1D::XSize)
      .def("Length",&Array1D::Length)
      .def("Resize",static_cast<void (Array1D::*)(const int&)>(&Array1D::Resize))
      .def("Resize",static_cast<void (Array1D::*)(const int&,const int&)>(&Array1D::Resize))
      .def("SetValue",&Array1D::SetValue)
      .def("PushBack",&Array1D::PushBack)
      .def("GetArrayPointer",&Array1D::GetArrayPointer)
      .def("GetConstArrayPointer",&Array1D::GetConstArrayPointer)
      .def("element",&Array1D::operator())
      .def("insert",static_cast<void (Array1D::*)(Array1D<int>&,int)>(&Array1D::insert))
      .def("insert",static_cast<void (Array1D::*)(const int&,int)>(&Array1D::insert))
      .def("erase",&Array1D::erase)
      .def("DumpBinary",static_cast<void (Array1D::*)(FILE*)>(&Array1D::DumpBinary)))
      .def("DumpBinary",static_cast<void (Array1D::*)(char*)>(&Array1D::DumpBinary)))
      .def("ReadBinary",static_cast<void (Array1D::*)(FILE*)>(&Array1D::ReadBinary)))
      .def("ReadBinary",static_cast<void (Array1D::*)(char*)>(&Array1D::ReadBinary)))
      .def("pyElement",&Array1D::operator[])
      .def("ReadBinary4py",&Array1D::ReadBinary4py)
      .def("DumpBinary4py",&Array1D::DumpBinary4py)
      .def("setArray",&Array1D::setArray)
      .def("setnpintArray",&Array1D::setnpintArray)
      .def("getnpintArray",&Array1D::getnpintArray)
      .def("flatten",&Array1D::flatten)
      .def("type",&Array1D::type)
      .def_property("xsize_",&Array1D::xsize_)
      .def_property("data_",&Array1D::data_)
      ;

      py::class_<Array1D<double>>(m, "Array1D<double>")
        .def(py::init<>())
        .def(py::init<const int&>())
        .def(py::init<const int&,const int&>())
        .def("Assign", &Array1D::operator=)
        .def(py::init<const Array1D &>())
        .def("Clear",&Array1D::Clear)
        .def("XSize",&Array1D::XSize)
        .def("Length",&Array1D::Length)
        .def("Resize",static_cast<void (Array1D::*)(const int&)>(&Array1D::Resize))
        .def("Resize",static_cast<void (Array1D::*)(const int&,const double&)>(&Array1D::Resize))
        .def("SetValue",&Array1D::SetValue)
        .def("PushBack",&Array1D::PushBack)
        .def("GetArrayPointer",&Array1D::GetArrayPointer)
        .def("GetConstArrayPointer",&Array1D::GetConstArrayPointer)
        .def("element",&Array1D::operator())
        .def("insert",static_cast<void (Array1D::*)(Array1D<double>&,int)>(&Array1D::insert))
        .def("insert",static_cast<void (Array1D::*)(const double&,int)>(&Array1D::insert))
        .def("erase",&Array1D::erase)
        .def("DumpBinary",static_cast<void (Array1D::*)(FILE*)>(&Array1D::DumpBinary)))
        .def("DumpBinary",static_cast<void (Array1D::*)(char*)>(&Array1D::DumpBinary)))
        .def("ReadBinary",static_cast<void (Array1D::*)(FILE*)>(&Array1D::ReadBinary)))
        .def("ReadBinary",static_cast<void (Array1D::*)(char*)>(&Array1D::ReadBinary)))
        .def("pyElement",&Array1D::operator[])
        .def("ReadBinary4py",&Array1D::ReadBinary4py)
        .def("DumpBinary4py",&Array1D::DumpBinary4py)
        .def("setArray",&Array1D::setArray)
        .def("setnpdblArray",&Array1D::setnpintArray)
        .def("getnpdblArray",&Array1D::getnpintArray)
        .def("flatten",&Array1D::flatten)
        .def("type",&Array1D::type)
        .def_property("xsize_",&Array1D::xsize_)
        .def_property("data_",&Array1D::data_)
        ;

      py::class_<Array2D<int>>(m,"Array2D")
        .def(py::init<>())
        .def(py::init<const int&,const int&>())
        .def(py::init<const int&,const int&>,const int&())
        .def(py::init<const Array2D &>())
        .def("Clear",&Array2D::Clear)
        .def("XSize",&Array2D::XSize)
        .def("YSize",&Array2D::YSize)
        .def("Resize",static_cast<void (Array2D::*)(const int&,const int&)>(&Array2D::Resize))
        .def("Resize",static_cast<void (Array2D::*)(const int&,const int&,const int&)>(&Array2D::Resize))
        .def("SetValue",&Array2D::SetValue)
        .def("GetArrayPointer",&Array2D::GetArrayPointer)
        .def("GetConstArrayPointer",&Array2D::GetConstArrayPointer)
        .def("element",&Array2D::operator())
        .def("insertRow",static_cast<void (Array2D::*)(Array1D<int>&,int)>(&Array2D::insertRow))
        .def("insertRow",static_cast<void (Array2D::*)(Array2D<int>&,int)>(&Array2D::insertRow))
        .def("eraseRow",&Array2D::eraseRow)
        .def("insertCol",static_cast<void (Array2D::*)(Array1D<int>&,int)>(&Array2D::insertCol))
        .def("insertCol",static_cast<void (Array2D::*)(Array2D<int>&,int)>(&Array2D::insertCol))
        .def("eraseCol",&Array2D::eraseCol)
        .def("DumpBinary",static_cast<void (Array2D::*)(FILE*)>(&Array2D::DumpBinary)))
        .def("DumpBinary",static_cast<void (Array2D::*)(char*)>(&Array2D::DumpBinary)))
        .def("ReadBinary",static_cast<void (Array2D::*)(FILE*)>(&Array2D::ReadBinary)))
        .def("ReadBinary",static_cast<void (Array2D::*)(char*)>(&Array2D::ReadBinary)))
        .def("pyElement",&Array1D::operator[])
        .def("getRow",&Array1D::getRow)
        .def("ReadBinary4py",&Array2D::ReadBinary4py)
        .def("DumpBinary4py",&Array2D::DumpBinary4py)
        .def("setArray",&Array2D::setArray)
        .def("flatten",&Array2D::flatten)
        .def("type",&Array2D::type)
        .def("setnpintArray",&Array2D::setnpintArray)
        .def("getnpintArray",&Array2D::getnpintArray)
        .def("setnpdblArray",&Array2D::setnpdblArray)
        .def("getnpdblArray",&Array2D::getnpdblArray)
        .def_property("xsize_",&Array2D::xsize_)
        .def_property("ysize_",&Array2D::xsize_)
        .def_property("data_",&Array2D::data_)
        .def_property("arraycopy",&Array2D::arraycopy)
        .def_property("rowvec",&Array2D::rowvec)
        ;

      py::class_<Array2D<double>>(m,"Array2D")
        .def(py::init<>())
        .def(py::init<const int&,const int&>())
        .def(py::init<const int&,const int&>,const int&())
        .def(py::init<const Array2D &>())
        .def("Clear",&Array2D::Clear)
        .def("XSize",&Array2D::XSize)
        .def("YSize",&Array2D::YSize)
        .def("Resize",static_cast<void (Array2D::*)(const int&,const int&)>(&Array2D::Resize))
        .def("Resize",static_cast<void (Array2D::*)(const int&,const int&,const double&)>(&Array2D::Resize))
        .def("SetValue",&Array2D::SetValue)
        .def("GetArrayPointer",&Array2D::GetArrayPointer)
        .def("GetConstArrayPointer",&Array2D::GetConstArrayPointer)
        .def("element",&Array2D::operator())
        .def("insertRow",static_cast<void (Array2D::*)(Array1D<double>&,int)>(&Array2D::insertRow))
        .def("insertRow",static_cast<void (Array2D::*)(Array2D<double>&,int)>(&Array2D::insertRow))
        .def("eraseRow",&Array2D::eraseRow)
        .def("insertCol",static_cast<void (Array2D::*)(Array1D<double>&,int)>(&Array2D::insertCol))
        .def("insertCol",static_cast<void (Array2D::*)(Array2D<double>&,int)>(&Array2D::insertCol))
        .def("eraseCol",&Array2D::eraseCol)
        .def("DumpBinary",static_cast<void (Array2D::*)(FILE*)>(&Array2D::DumpBinary)))
        .def("DumpBinary",static_cast<void (Array2D::*)(char*)>(&Array2D::DumpBinary)))
        .def("ReadBinary",static_cast<void (Array2D::*)(FILE*)>(&Array2D::ReadBinary)))
        .def("ReadBinary",static_cast<void (Array2D::*)(char*)>(&Array2D::ReadBinary)))
        .def("pyElement",&Array1D::operator[])
        .def("getRow",&Array1D::getRow)
        .def("ReadBinary4py",&Array2D::ReadBinary4py)
        .def("DumpBinary4py",&Array2D::DumpBinary4py)
        .def("setArray",&Array2D::setArray)
        .def("flatten",&Array2D::flatten)
        .def("type",&Array2D::type)
        .def("setnpintArray",&Array2D::setnpintArray)
        .def("getnpintArray",&Array2D::getnpintArray)
        .def("setnpdblArray",&Array2D::setnpdblArray)
        .def("getnpdblArray",&Array2D::getnpdblArray)
        .def_property("xsize_",&Array2D::xsize_)
        .def_property("ysize_",&Array2D::xsize_)
        .def_property("data_",&Array2D::data_)
        .def_property("arraycopy",&Array2D::arraycopy)
        .def_property("rowvec",&Array2D::rowvec)
        ;

      m.def("read_datafile",&read_datafile)
      m.def("read_datafileVS",&read_datafileVS)
      m.def("read_datafile_1d",&read_datafile_1d)
      m.def("write_datafile_size",&write_datafile_size)
      m.def("write_datafile",&write_datafile)
      m.def("write_datafile_1d",&write_datafile_1d)

      m.def("array1Dto2D",&array1Dto2D)
      m.def("array2Dto1D",&array2Dto1D)
      m.def("paste",&paste)
      m.def("generate_multigrid",&generate_multigrid)
      m.def("merge",&merge)
      m.def("append",&append)
      m.def("transpose",&transpose)
      m.def("flatten",&flatten)
      m.def("fold_1dto2d_rowfirst",&fold_1dto2d_rowfirst)
      m.def("fold_1dto2d_colfirst",&fold_1dto2d_colfirst)
      m.def("swap",&swap)
      m.def("access",&access)
      m.def("getRow",&getRow)
      m.def("getCol",&getCol)
      m.def("addVal",&addVal)
      m.def("subVector",&subVector)
      m.def("subMatrix_row",&subMatrix_row)
      m.def("subMatrix_col",&subMatrix_col)
      m.def("matPvec",&matPvec)
      m.def("maxVal",&maxVal)
      m.def("setdiff",&setdiff)
      m.def("setdiff_s",&setdiff_s)
      m.def("shell_sort",&shell_sort)
      m.def("shell_sort_col",&shell_sort_col)
      m.def("shell_sort_all",&shell_sort_all)
      m.def("quicksort3",&quicksort3)
      m.def("intersect",&intersect)
      m.def("find",&find)
      m.def("prodAlphaMatVec",&prodAlphaMatVec)
      m.def("prodAlphaMatTVec",&prodAlphaMatTVec)
      m.def("prodAlphaMatMat",&prodAlphaMatMat)
      m.def("prodAlphaMatTMat",&prodAlphaMatTMat)
      m.def("addVecAlphaVecPow",&addVecAlphaVecPow)
      m.def("prod_vecTmatvec",&prod_vecTmatvec)
      m.def("MatTMat",&MatTMat)
      m.def("delRow",&delRow)
      m.def("delCol",&delCol)
      m.def("paddMatRow",&paddMatRow)
      m.def("paddMatCol",&paddMatCol)
      m.def("paddMatRow",&paddMatColScal)
      m.def("is_equal",&is_equal)
      m.def("is_less",&is_less)
      m.def("vecIsInArray",&vecIsInArray)
      m.def("select_kth",&select_kth)
      m.def("logdeterm",&logdeterm)
      m.def("trace",&trace)
      m.def("evalLogMVN",&evalLogMVN)
      m.def("diag",&diag)
      m.def("copy",&copy)
      m.def("mtxdel",&mtxdel)
      m.def("add",&add)
      m.def("addinplace",&addinplace)
      m.def("subtract",&subtract)
      m.def("subtractinplace",&subtractinplace)
      m.def("scale",&scale)
      m.def("scaleinplace",&scaleinplace)
      m.def("dotdivide",&dotdivide)
      m.def("dotmult",&dotmult)
      m.def("norm",&norm)
      m.def("dist_sq",&dist_sq)
      m.def("Trans",&Trans)
      m.def("dot",&dot)
      m.def("dotT",&dotT)
      m.def("INV",&INV)
      m.def("AinvH",&AinvH)
      m.def("Ainvb",&Ainvb)
      m.def("LSTSQ",&LSTSQ)
      m.def("QR",&QR)
      m.def("SVD",&SVD)
      m.def("printarray",printarray)





}
