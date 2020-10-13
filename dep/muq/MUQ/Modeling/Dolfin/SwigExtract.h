#ifndef SWIGEXTRACT_H_
#define SWIGEXTRACT_H_

#include <pybind11/pybind11.h>

struct PySwigObject {
    PyObject_HEAD
    void * ptr;
    const char * desc;
};


namespace muq{
namespace Modeling{

/** @class SwigExtract
    @ingroup FenicsCoupling

    @brief Class to extract a c++ class that was exposed to python with Swig.
    @detials This class is used to unwrap a c++ class that was exposed to python with Swig and then passed to MUQ through pybind11 bindings.
             The class is defined in a way to enable straightforward conversions with code like
@code
  pybind11::object pyObject = ... ;
  std::shared_ptr<CppClass> cppPtr = SwigExtract(pyObject);
@endcode

Instead of using the implicit conversion, it is also possible to explicitly cast the python object:
@code
  pybind11::object pyObject = ... ;
  SwigExtract extractor(pyObject);
  std::shared_ptr<CppClass> cppPtr = extractor.Cast<std::shared_pointer<CppClass>>();
@endcode 

*/
    class SwigExtract
    {
        
    public:
        SwigExtract(pybind11::object const& objIn) : obj(objIn){};
        
        template<typename T>
        operator T()
        {
            return Cast<T>();
        }
        
        template<typename T>
        T Cast()
        {
            // This functions performs some magic performed by the boost.python folks and described here:
            // https://wiki.python.org/moin/boost.python/HowTo#SWIG_exposed_C.2B-.2B-_object_from_Python
            
            PyObject* objPtr = pybind11::handle(obj).ptr();
            
            char thisStr[] = "this";
            
            //first we need to get the this attribute from the Python Object
            if (!PyObject_HasAttrString(objPtr, thisStr))
                ConversionError<T>();
            
            PyObject* thisAttr = PyObject_GetAttrString(objPtr, thisStr);
            if (thisAttr == NULL)
                ConversionError<T>();   
            
            T* pointer = (T*) ((PySwigObject*)thisAttr)->ptr;
            Py_DECREF(thisAttr);
            
            return *pointer;
        }
        
    private:
        
        pybind11::object const& obj;
        
        template<typename T>
        void ConversionError()
        {   
            throw std::invalid_argument("Unable to convert Python type of " + std::string(pybind11::str(obj.attr("__repr__"))) + " to c++ type of " + std::string(typeid(T).name()));
        }
    };


} // namepsace Modeling
} // namespace muq

#endif
