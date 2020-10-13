#ifndef BOOSTANYSERIALIZER_H
#define BOOSTANYSERIALIZER_H

#include "cereal/cereal.hpp"

#include "parcer/Eigen.h"
#include <boost/any.hpp>
#include <cereal/types/string.hpp>

#include <typeinfo>
#include <string>

namespace cereal{

struct TemporaryBoostAnyConst
{
  TemporaryBoostAnyConst(boost::any const& anyIn) : any(anyIn) {};
  boost::any const& any;
};
struct TemporaryBoostAny
{
  TemporaryBoostAny() = default;
  TemporaryBoostAny(boost::any& anyIn) : any(anyIn) {};
  boost::any any;
};

template<class Archive>
class BoostAnySerializer
{
public:

  // We assume that this function is called after the type name has been read.
  // Thus, only the object has to be loaded from the archive.
  template<typename BaseType>
  static boost::any load(Archive & ar)
  {
    BaseType temp;
    ar(temp);
    return boost::any(temp);
  };

  template<typename BaseType>
  static void save(boost::any obj, Archive & ar)
  {
    // First, save the type name
    const std::string& typeName = obj.type().name();
    ar(typeName);

    // Now, save the object
    BaseType const& val = boost::any_cast<const BaseType&>(obj);
    ar(val);
  };
};

/**
Serializes a boost::any object.  Assumes that the boost::any wraps around a double, Eigen::VectorXd, or Eigen::MatrixXd
*/
template<class Archive>
void save(Archive & ar, TemporaryBoostAnyConst const& obj)
{
  if(obj.any.type() == typeid(double)){
    BoostAnySerializer<Archive>::template save<double>(obj.any,ar);
  }else if(obj.any.type() == typeid(float)){
    BoostAnySerializer<Archive>::template save<float>(obj.any,ar);
  }else if(obj.any.type() == typeid(std::string)){
    BoostAnySerializer<Archive>::template save<std::string>(obj.any,ar);
  }else if(obj.any.type() == typeid(int)){
    BoostAnySerializer<Archive>::template save<int>(obj.any,ar);
  }else if(obj.any.type() == typeid(unsigned int)){
    BoostAnySerializer<Archive>::template save<unsigned int>(obj.any,ar);
    
  }else if(obj.any.type() == typeid(Eigen::Vector2d)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector2d>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Vector3d)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector3d>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Vector4d)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector4d>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::VectorXd)){
    BoostAnySerializer<Archive>::template save<Eigen::VectorXd>(obj.any,ar);
    
  }else if(obj.any.type() == typeid(Eigen::Vector2f)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector2f>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Vector3f)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector3f>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Vector4f)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector4f>(obj.any,ar);
      }else if(obj.any.type() == typeid(Eigen::VectorXf)){
    BoostAnySerializer<Archive>::template save<Eigen::VectorXf>(obj.any,ar);
    
  }else if(obj.any.type() == typeid(Eigen::Vector2i)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector2i>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Vector3i)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector3i>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Vector4i)){
    BoostAnySerializer<Archive>::template save<Eigen::Vector4i>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::VectorXi)){
    BoostAnySerializer<Archive>::template save<Eigen::VectorXi>(obj.any,ar);
    
  }else if(obj.any.type() == typeid(Eigen::Matrix2d)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix2d>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Matrix3d)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix3d>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Matrix4d)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix4d>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::MatrixXd)){
    BoostAnySerializer<Archive>::template save<Eigen::MatrixXd>(obj.any,ar);
    
  }else if(obj.any.type() == typeid(Eigen::Matrix2f)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix2f>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Matrix3f)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix3f>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Matrix4f)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix4f>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::MatrixXf)){
    BoostAnySerializer<Archive>::template save<Eigen::MatrixXf>(obj.any,ar);
    
  }else if(obj.any.type() == typeid(Eigen::Matrix2i)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix2i>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Matrix3i)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix3i>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::Matrix4i)){
    BoostAnySerializer<Archive>::template save<Eigen::Matrix4i>(obj.any,ar);
  }else if(obj.any.type() == typeid(Eigen::MatrixXi)){
    BoostAnySerializer<Archive>::template save<Eigen::MatrixXi>(obj.any,ar);
  }else{
    assert(false);
  }
};

/**
  Deserializes a boost any.
*/
template<class Archive>
void load(Archive & ar, boost::any & obj)
{
  // First, read the type
  std::string typeName;
  ar(typeName);

  // Now read the value itself and store in a boost::any
  if(typeName == typeid(double).name()){
    obj = BoostAnySerializer<Archive>::template load<double>(ar);
  }else if(typeName == typeid(float).name()){
    obj = BoostAnySerializer<Archive>::template load<float>(ar);
  }else if(typeName == typeid(std::string).name()){
    obj = BoostAnySerializer<Archive>::template load<std::string>(ar);
  }else if(typeName == typeid(int).name()){
    obj = BoostAnySerializer<Archive>::template load<int>(ar);
  }else if(typeName == typeid(unsigned int).name()){
    obj = BoostAnySerializer<Archive>::template load<unsigned int>(ar);

  }else if(typeName == typeid(Eigen::Vector2d).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector2d>(ar);
  }else if(typeName == typeid(Eigen::Vector3d).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector3d>(ar);
  }else if(typeName == typeid(Eigen::Vector4d).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector4d>(ar);
  }else if(typeName == typeid(Eigen::VectorXd).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::VectorXd>(ar);
    
  }else if(typeName == typeid(Eigen::Vector2f).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector2f>(ar);
  }else if(typeName == typeid(Eigen::Vector3f).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector3f>(ar);
  }else if(typeName == typeid(Eigen::Vector4f).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector4f>(ar);
  }else if(typeName == typeid(Eigen::VectorXf).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::VectorXf>(ar);
    
  }else if(typeName == typeid(Eigen::Vector2i).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector2i>(ar);
  }else if(typeName == typeid(Eigen::Vector3i).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector3i>(ar);
  }else if(typeName == typeid(Eigen::Vector4i).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Vector4i>(ar);
  }else if(typeName == typeid(Eigen::VectorXi).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::VectorXi>(ar);
    
  }else if(typeName == typeid(Eigen::Matrix2d).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix2d>(ar);
  }else if(typeName == typeid(Eigen::Matrix3d).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix3d>(ar);
  }else if(typeName == typeid(Eigen::Matrix4d).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix4d>(ar);
  }else if(typeName == typeid(Eigen::MatrixXd).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::MatrixXd>(ar);
    
  }else if(typeName == typeid(Eigen::Matrix2f).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix2f>(ar);
  }else if(typeName == typeid(Eigen::Matrix3f).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix3f>(ar);
  }else if(typeName == typeid(Eigen::Matrix4f).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix4f>(ar);
  }else if(typeName == typeid(Eigen::MatrixXf).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::MatrixXf>(ar);
    
  }else if(typeName == typeid(Eigen::Matrix2i).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix2i>(ar);
  }else if(typeName == typeid(Eigen::Matrix3i).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix3i>(ar);
  }else if(typeName == typeid(Eigen::Matrix4i).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::Matrix4i>(ar);
  }else if(typeName == typeid(Eigen::MatrixXi).name()){
    obj = BoostAnySerializer<Archive>::template load<Eigen::MatrixXi>(ar);
  }else{
    assert(false);
  }
}

} // namespace cereal

#endif
