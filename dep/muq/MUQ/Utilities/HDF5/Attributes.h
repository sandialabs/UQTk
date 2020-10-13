#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

#include "MUQ/Utilities/HDF5/HDF5File.h"

#include <map>
#include <string>

namespace muq
{

namespace Utilities
{
    
    
class Attribute
{

public:

    Attribute(){};
    
    Attribute(std::shared_ptr<HDF5File> file_,
	      std::string const& path_,
	      std::string const& name_) : file(file_),
	                                  path(path_),
	                                  name(name_){};

    virtual ~Attribute(){};
  
  
  //////////////////////////////////////////////////////
  // conversions FROM other types
  //////////////////////////////////////////////////////
  template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
  Attribute& operator=(ScalarType val)
  {
      assert(file);
      file->WriteScalarAttribute(path, name, val);
      return *this;
  }

  template<typename ScalarType, int fixedRows>
  Attribute& operator=(Eigen::Matrix<ScalarType, fixedRows, 1> const& val)
  {
      assert(file);
      file->WriteVectorAttribute(path, name, val);
      return *this;
  }
  
  Attribute& operator=(std::string const& val);
  
  //////////////////////////////////////////////////////
  // conversions TO other types
  //////////////////////////////////////////////////////
  template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
  operator ScalarType() const
  {
      assert(file);
       
      return file->GetScalarAttribute<ScalarType>(path,name);
  };

  template<typename ScalarType>
  operator Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>() const
  {
      assert(file);
      return file->GetVectorAttribute<ScalarType>(path,name);
  };
  
  operator std::string() const;
  
private:

  std::shared_ptr<HDF5File> file;
  
  std::string path; // the path into the HDF5 file
  std::string name; // the name of the attribute
  
};

class AttributeList
{

public:

    AttributeList() {};
    
    AttributeList(std::shared_ptr<HDF5File>        file_,
		  std::string               const& path_) : file(file_),
                                                            path(path_){};
    
    Attribute& operator[](std::string const& attrName);

    
private:
    std::map<std::string, Attribute> attributes;
    
    std::shared_ptr<HDF5File> file;
    
    std::string path;
};

} // namespace Utilities
} // namespace muq

#endif
