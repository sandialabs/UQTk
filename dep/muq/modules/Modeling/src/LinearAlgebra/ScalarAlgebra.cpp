#include "MUQ/Modeling/LinearAlgebra/ScalarAlgebra.h"

using namespace muq::Modeling;

ScalarAlgebra::ScalarAlgebra() {}

ScalarAlgebra::~ScalarAlgebra() {}

bool ScalarAlgebra::IsScalar(std::type_info const& obj_type) {
  // is this a scalar type?
  return typeid(double)==obj_type || typeid(float)==obj_type || typeid(int)==obj_type || typeid(unsigned int)==obj_type;
}

bool ScalarAlgebra::IsZero(boost::any const& obj) {
  if( obj.type()==typeid(double) ) { return boost::any_cast<double const>(obj)==0.0; }
  if( obj.type()==typeid(float) ) { return boost::any_cast<float const>(obj)==0.0; }
  if( obj.type()==typeid(int) ) { return boost::any_cast<int const>(obj)==0; }
  if( obj.type()==typeid(unsigned int) ) { return boost::any_cast<unsigned int const>(obj)==0; }

  // something when wrong
  assert(false);
  return false;
}

boost::any ScalarAlgebra::Zero(std::type_info const& type) {
  if( typeid(double)==type ) { return (double)0.0; }
  if( typeid(float)==type ) { return (float)0.0; }
  if( typeid(int)==type ) { return (int)0; }
  if( typeid(unsigned int)==type ) { return (unsigned int)0; }

  // something went wrong
  assert(false);
  return boost::none;
}

double ScalarAlgebra::Norm(boost::any const& obj) {
  if( typeid(double)==obj.type() ) { return Magnitude<double>(obj); }
  if( typeid(float)==obj.type() ) { return Magnitude<float>(obj); }
  if( typeid(int)==obj.type() ) { return Magnitude<int>(obj); }
  if( typeid(unsigned int)==obj.type() ) { return Magnitude<unsigned int>(obj); }

  // something went wrong
  assert(false);
  return -1.0;
}

double ScalarAlgebra::InnerProduct(boost::any const& in0, boost::any const& in1) {
  const boost::any result = ScalarAlgebra::Multiply(in0, in1);

  if( result.type()==typeid(double) ) { return boost::any_cast<double const>(result); }
  if( result.type()==typeid(float) ) { return (double)boost::any_cast<float const>(result); }
  if( result.type()==typeid(int) ) { return (double)boost::any_cast<int const>(result); }
  if( result.type()==typeid(unsigned int) ) { return (double)boost::any_cast<unsigned int const>(result); }

  // something went wrong
  return std::numeric_limits<double>::quiet_NaN();
}

boost::any ScalarAlgebra::OuterProduct(boost::any const& in0, boost::any const& in1) {
  return ScalarAlgebra::Multiply(in0, in1);
  /*if( in0.type()==typeid(double) ) { return Multiply<double>(in0, in1); }
  if( in0.type()==typeid(float) ) { return Multiply<float>(in0, in1); }
  if( in0.type()==typeid(int) ) { return Multiply<int>(in0, in1); }
  if( in0.type()==typeid(unsigned int) ) { return Multiply<unsigned int>(in0, in1); }

  // something went wrong
  assert(false);
  return boost::none;*/
}

boost::any ScalarAlgebra::Identity(std::type_info const& type) {
  if( type==typeid(double) ) { return (double)1.0; }
  if( type==typeid(float) ) { return (float)1.0; }
  if( type==typeid(int) ) { return (int)1; }
  if( type==typeid(unsigned int) ) { return (unsigned int)1; }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any ScalarAlgebra::Add(boost::any const& in0, boost::any const& in1) {
  if( in0.type()==typeid(double) ) { return Add<double>(in0, in1); }
  if( in0.type()==typeid(float) ) { return Add<float>(in0, in1); }
  if( in0.type()==typeid(int) ) { return Add<int>(in0, in1); }
  if( in0.type()==typeid(unsigned int) ) { return Add<unsigned int>(in0, in1); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any ScalarAlgebra::Subtract(boost::any const& in0, boost::any const& in1) {
  if( in0.type()==typeid(double) ) { return Subtract<double>(in0, in1); }
  if( in0.type()==typeid(float) ) { return Subtract<float>(in0, in1); }
  if( in0.type()==typeid(int) ) { return Subtract<int>(in0, in1); }
  if( in0.type()==typeid(unsigned int) ) { return Subtract<unsigned int>(in0, in1); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any ScalarAlgebra::Inverse(boost::any const& obj) {
  if( obj.type()==typeid(double) ) { return 1.0/boost::any_cast<double const>(obj); }
  if( obj.type()==typeid(float) ) { return (float)1.0/boost::any_cast<float const>(obj); }
  if( obj.type()==typeid(int) ) { return 1.0/boost::any_cast<int const>(obj); }
  if( obj.type()==typeid(unsigned int) ) { return 1.0/boost::any_cast<unsigned int const>(obj); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any ScalarAlgebra::Multiply(boost::any const& in0, boost::any const& in1) {
  if( in0.type()==typeid(double) ) { return Multiply<double>(in0, in1); }
  if( in0.type()==typeid(float) ) { return Multiply<float>(in0, in1); }
  if( in0.type()==typeid(int) ) { return Multiply<int>(in0, in1); }
  if( in0.type()==typeid(unsigned int) ) { return Multiply<unsigned int>(in0, in1); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any ScalarAlgebra::SquareRoot(boost::any const& obj) {
  if( obj.type()==typeid(double) ) { return std::sqrt(boost::any_cast<double>(obj)); }
  if( obj.type()==typeid(float) ) { return std::sqrt(boost::any_cast<float>(obj)); }
  if( obj.type()==typeid(int) ) { return std::sqrt(boost::any_cast<int>(obj)); }
  if( obj.type()==typeid(unsigned int) ) { return std::sqrt(boost::any_cast<unsigned int>(obj)); }

  // something went wrong
  assert(false);
  return boost::none;
}

double ScalarAlgebra::LogDeterminate(boost::any const& obj) {
  if( typeid(double)==obj.type() ) { return std::log(Magnitude<double>(obj)); }
  if( typeid(float)==obj.type() ) { return std::log(Magnitude<float>(obj)); }
  if( typeid(int)==obj.type() ) { return std::log(Magnitude<int>(obj)); }
  if( typeid(unsigned int)==obj.type() ) { return std::log(Magnitude<unsigned int>(obj)); }

  // something went wrong
  assert(false);
  return -1.0;
}
