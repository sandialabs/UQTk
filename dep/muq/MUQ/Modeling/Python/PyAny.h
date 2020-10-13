#ifndef PYANY_H_
#define PYANY_H_

#include <boost/any.hpp>

#include <typeindex>

#include <cxxabi.h>

#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"

namespace pybind11 { namespace detail {

    template <> struct type_caster<boost::any> {
    public:
        
        PYBIND11_TYPE_CASTER(boost::any, _("boost::any"));

        bool load(handle src, bool) {

            PyObject *source = src.ptr();

            char thisStr[] = "this";

            // If the pyobject does not have a "this" attribute, then we should try to convert it to a native c++ type (e.g., double, int, string, etc...)
            if (!PyObject_HasAttrString(source, thisStr))
            {
                std::string typeStr =  source->ob_type->tp_name;
                auto it = fromPythonMap.find(typeStr);
                if( it != fromPythonMap.end()){
                    value = it->second(source);
                    return true;
                }else{
                    value = boost::any(src.cast<pybind11::object>());
                    return true;
                }
            }

            // In this case, the python object is a class
            PyObject* thisAttr = PyObject_GetAttrString(source, thisStr);
            if (thisAttr == NULL)
            {
                std::cerr << "\nERROR: Could not get the this attr.\n\n";
            }

            std::cerr << "ERROR: Don't have a good way to deal with generic classes yet...." << std::endl;
            assert(false);
            
            Py_DECREF(thisAttr);
            return false;
        }

        static handle cast(boost::any src, return_value_policy /* policy */, handle /* parent */) {

            auto it = toPythonMap.find(src.type());
            if(it != toPythonMap.end()){
                return it->second(src);
            }else{
                std::cerr << "WARNING: Could not convert type " << demangle_typename(src.type().name()) << " directly to Python type.  Trying to cast as pybind11::object." << std::endl;
                return boost::any_cast<pybind11::object>(src);
            }
        }

    private:

        using EigenRowType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using EigenColType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

        using EigenVectorType = Eigen::Matrix<double, Eigen::Dynamic,1>;
        using EigenRowVectorType = Eigen::Matrix<double, 1, Eigen::Dynamic>;
        
        // A map of functions to convert a boost::any to a corresponding PyObject
        static std::map<std::type_index, std::function<handle(boost::any const&)>> toPythonMap;
        static std::map<std::string, std::function<boost::any(PyObject*)>> fromPythonMap;
        
        static std::map<std::type_index, std::function<handle(boost::any const&)>> createToPythonMap()
        {
          std::map<std::type_index, std::function<handle(boost::any const&)>> m;
          m[typeid(long)]          = [](boost::any const& x) { return PyLong_FromLong(boost::any_cast<long>(x));};
          m[typeid(int)]           = [](boost::any const& x) { return PyLong_FromLong(boost::any_cast<int>(x));};
          m[typeid(unsigned)]      = [](boost::any const& x) { return PyLong_FromLong(boost::any_cast<unsigned>(x));};
          m[typeid(unsigned long)] = [](boost::any const& x) { return PyLong_FromLong(boost::any_cast<unsigned long>(x));};
          m[typeid(double)]        = [](boost::any const& x) { return PyFloat_FromDouble(boost::any_cast<double>(x));};
          m[typeid(float)]         = [](boost::any const& x) { return PyFloat_FromDouble(static_cast<double>(boost::any_cast<float>(x)));};
          m[typeid(std::string)]   = [](boost::any const& x) { return PyBytes_FromString(boost::any_cast<std::string>(x).c_str());};
          
          m[typeid(Eigen::Ref<EigenRowType>)] = [](boost::any const& x){
              return eigen_array_cast<EigenProps<Eigen::Ref<EigenRowType>>>(boost::any_cast<Eigen::Ref<EigenRowType>>(x));
          };

          m[typeid(Eigen::Ref<EigenColType>)] = [](boost::any const& x){
              return eigen_array_cast<EigenProps<Eigen::Ref<EigenColType>>>(boost::any_cast<Eigen::Ref<EigenColType>>(x));
          };

          m[typeid(EigenRowType)] = [](boost::any const& x){
              return eigen_array_cast<EigenProps<EigenRowType>>(boost::any_cast<EigenRowType>(x));
          };

          m[typeid(EigenColType)] = [](boost::any const& x){
              return eigen_array_cast<EigenProps<EigenColType>>(boost::any_cast<EigenColType>(x));
          };
          
          m[typeid(EigenVectorType)] = [](boost::any const& x){
              return eigen_array_cast<EigenProps<EigenVectorType>>(boost::any_cast<EigenVectorType>(x));
          };

          m[typeid(EigenRowVectorType)] = [](boost::any const& x){
              return eigen_array_cast<EigenProps<EigenRowVectorType>>(boost::any_cast<EigenRowVectorType>(x));
          };

                    
          return m;
        }

        static std::map<std::string, std::function<boost::any(PyObject*)>> createFromPythonMap()
        {
          std::map<std::string, std::function<boost::any(PyObject*)>> m;
          m["float"] = [](PyObject* obj) { return boost::any(PyFloat_AsDouble(obj));};
          m["int"] = [](PyObject* obj) { return boost::any(PyLong_AsLong(obj));};
          m["str"] = [](PyObject* obj){
              PyObject* str_exc_type = PyObject_Repr(obj);
                
              PyObject* pyStr = PyUnicode_AsEncodedString(str_exc_type, "utf-8", "Error ~");
              const char *strExcType =  PyBytes_AS_STRING(pyStr);
              
              boost::any output = std::string(strExcType);
              
              Py_XDECREF(str_exc_type);
              Py_XDECREF(pyStr);
              
              return output;
          };
          
          m["numpy.ndarray"] = [](PyObject* obj){
              
              type_caster<Eigen::Ref<EigenRowType>> rowCaster;
              bool res = rowCaster.load(obj, false);
              if(res){
                  Eigen::Ref<EigenRowType> *ref = rowCaster;
                  return boost::any(*ref);
              }else{
                  type_caster<Eigen::Ref<EigenColType>> colCaster;
                  bool res = colCaster.load(obj, false);

                  if(!res){
                      std::cerr << "ERROR: Could not convert numpy array to Eigen::Ref.  Current support is only for row-major and col-major arrays of doubles.  Is the numpy array full of doubles?";
                      assert(res);
                  }

                  Eigen::Ref<EigenColType> *ref = colCaster;
                  return boost::any(*ref);
              }
          };
          
          return m;
        }



        static std::string demangle_typename(const char* name) {
            
            int status = -4; // some arbitrary value to eliminate the compiler warning
            
            // enable c++11 by passing the flag -std=c++11 to g++
            std::unique_ptr<char, void(*)(void*)> res {
                abi::__cxa_demangle(name, NULL, NULL, &status),
                    std::free
                    };
            
            return (status==0) ? res.get() : name ;
        }
    };

    // Fill in the map for functions that transform c++ types to python types
    std::map<std::type_index, std::function<handle(boost::any const&)>>
    type_caster<boost::any>::toPythonMap = type_caster<boost::any>::createToPythonMap();

    // Fill in the map for functions that transform python types to c++ types
    std::map<std::string, std::function<boost::any(PyObject*)>>
    type_caster<boost::any>::fromPythonMap = type_caster<boost::any>::createFromPythonMap();

}} // namespace pybind11::detail

#endif
