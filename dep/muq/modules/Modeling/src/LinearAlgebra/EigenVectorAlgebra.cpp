#include "MUQ/Modeling/LinearAlgebra/EigenVectorAlgebra.h"

using namespace muq::Modeling;

EigenVectorAlgebra::EigenVectorAlgebra() {}

EigenVectorAlgebra::~EigenVectorAlgebra() {}

bool EigenVectorAlgebra::IsEigenVector(std::type_info const& obj_type) {
  // is this an Eigen::Vector type?
  return typeid(Eigen::Vector2d)==obj_type
    || typeid(Eigen::Vector2f)==obj_type
    || typeid(Eigen::Vector2i)==obj_type
    || typeid(Eigen::Vector3d)==obj_type
    || typeid(Eigen::Vector3f)==obj_type
    || typeid(Eigen::Vector3i)==obj_type
    || typeid(Eigen::Vector4d)==obj_type
    || typeid(Eigen::Vector4f)==obj_type
    || typeid(Eigen::Vector4i)==obj_type
    || typeid(Eigen::VectorXd)==obj_type
    || typeid(Eigen::VectorXf)==obj_type
    || typeid(Eigen::VectorXi)==obj_type;
}

bool EigenVectorAlgebra::IsZero(boost::any const& obj) {
  if( typeid(Eigen::Vector2d)==obj.type() ) { return IsZero<Eigen::Vector2d>(obj); }
  if( typeid(Eigen::Vector2f)==obj.type() ) { return IsZero<Eigen::Vector2f>(obj); }
  if( typeid(Eigen::Vector2i)==obj.type() ) { return IsZero<Eigen::Vector2i>(obj); }

  if( typeid(Eigen::Vector3d)==obj.type() ) { return IsZero<Eigen::Vector3d>(obj); }
  if( typeid(Eigen::Vector3f)==obj.type() ) { return IsZero<Eigen::Vector3f>(obj); }
  if( typeid(Eigen::Vector3i)==obj.type() ) { return IsZero<Eigen::Vector3i>(obj); }

  if( typeid(Eigen::Vector4d)==obj.type() ) { return IsZero<Eigen::Vector4d>(obj); }
  if( typeid(Eigen::Vector4f)==obj.type() ) { return IsZero<Eigen::Vector4f>(obj); }
  if( typeid(Eigen::Vector4i)==obj.type() ) { return IsZero<Eigen::Vector4i>(obj); }

  if( typeid(Eigen::VectorXd)==obj.type() ) { return IsZero<Eigen::VectorXd>(obj); }
  if( typeid(Eigen::VectorXf)==obj.type() ) { return IsZero<Eigen::VectorXf>(obj); }
  if( typeid(Eigen::VectorXi)==obj.type() ) { return IsZero<Eigen::VectorXi>(obj); }

  // something went wrong
  assert(false);
  return false;
}

unsigned int EigenVectorAlgebra::Size(boost::any const& vec) {
  if( typeid(Eigen::Vector2d)==vec.type() ) { return 2; }
  if( typeid(Eigen::Vector3d)==vec.type() ) { return 3; }
  if( typeid(Eigen::Vector4d)==vec.type() ) { return 4; }
  if( typeid(Eigen::VectorXd)==vec.type() ) {
    const Eigen::VectorXd& eig = boost::any_cast<Eigen::VectorXd const&>(vec);
    return eig.size();
  }

  if( typeid(Eigen::Vector2f)==vec.type() ) { return 2; }
  if( typeid(Eigen::Vector3f)==vec.type() ) { return 3; }
  if( typeid(Eigen::Vector4f)==vec.type() ) { return 4; }
  if( typeid(Eigen::VectorXf)==vec.type() ) {
    const Eigen::VectorXf& eig = boost::any_cast<Eigen::VectorXf const&>(vec);
    return eig.size();
  }

  if( typeid(Eigen::Vector2i)==vec.type() ) { return 2; }
  if( typeid(Eigen::Vector3i)==vec.type() ) { return 3; }
  if( typeid(Eigen::Vector4i)==vec.type() ) { return 4; }
  if( typeid(Eigen::VectorXi)==vec.type() ) {
    const Eigen::VectorXi& eig = boost::any_cast<Eigen::VectorXi const&>(vec);
    return eig.size();
  }

  // something went wront
  assert(false);
  return 0;
}

double EigenVectorAlgebra::Norm(boost::any const& vec) {
  if( typeid(Eigen::Vector2d)==vec.type() ) { return boost::any_cast<Eigen::Vector2d const&>(vec).norm(); }
  if( typeid(Eigen::Vector3d)==vec.type() ) { return boost::any_cast<Eigen::Vector3d const&>(vec).norm(); }
  if( typeid(Eigen::Vector4d)==vec.type() ) { return boost::any_cast<Eigen::Vector4d const&>(vec).norm(); }
  if( typeid(Eigen::VectorXd)==vec.type() ) { return boost::any_cast<Eigen::VectorXd const&>(vec).norm(); }

  if( typeid(Eigen::Vector2f)==vec.type() ) { return boost::any_cast<Eigen::Vector2f const&>(vec).norm(); }
  if( typeid(Eigen::Vector3f)==vec.type() ) { return boost::any_cast<Eigen::Vector3f const&>(vec).norm(); }
  if( typeid(Eigen::Vector4f)==vec.type() ) { return boost::any_cast<Eigen::Vector4f const&>(vec).norm(); }
  if( typeid(Eigen::VectorXf)==vec.type() ) { return boost::any_cast<Eigen::VectorXf const&>(vec).norm(); }

  if( typeid(Eigen::Vector2i)==vec.type() ) { return boost::any_cast<Eigen::Vector2i const&>(vec).norm(); }
  if( typeid(Eigen::Vector3i)==vec.type() ) { return boost::any_cast<Eigen::Vector3i const&>(vec).norm(); }
  if( typeid(Eigen::Vector4i)==vec.type() ) { return boost::any_cast<Eigen::Vector4i const&>(vec).norm(); }
  if( typeid(Eigen::VectorXi)==vec.type() ) { return boost::any_cast<Eigen::VectorXi const&>(vec).norm(); }

  // something went wront
  assert(false);
  return 0;
}

double EigenVectorAlgebra::InnerProduct(boost::any const& vec1, boost::any const& vec2) {
  if( typeid(Eigen::Vector2d)==vec1.type() && typeid(Eigen::Vector2d)==vec2.type() ) { return InnerProduct<Eigen::Vector2d, Eigen::Vector2d>(vec1, vec2); }
  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::Vector2d)==vec2.type() ) { return InnerProduct<Eigen::VectorXd, Eigen::Vector2d>(vec1, vec2); }
  if( typeid(Eigen::Vector2d)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return InnerProduct<Eigen::Vector2d, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::Vector3d)==vec1.type() && typeid(Eigen::Vector3d)==vec2.type() ) { return InnerProduct<Eigen::Vector3d, Eigen::Vector3d>(vec1, vec2); }
  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::Vector3d)==vec2.type() ) { return InnerProduct<Eigen::VectorXd, Eigen::Vector3d>(vec1, vec2); }
  if( typeid(Eigen::Vector3d)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return InnerProduct<Eigen::Vector3d, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::Vector4d)==vec1.type() && typeid(Eigen::Vector4d)==vec2.type() ) { return InnerProduct<Eigen::Vector4d, Eigen::Vector4d>(vec1, vec2); }
  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::Vector4d)==vec2.type() ) { return InnerProduct<Eigen::VectorXd, Eigen::Vector4d>(vec1, vec2); }
  if( typeid(Eigen::Vector4d)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return InnerProduct<Eigen::Vector4d, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return InnerProduct<Eigen::VectorXd, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::Vector2f)==vec1.type() && typeid(Eigen::Vector2f)==vec2.type() ) { return InnerProduct<Eigen::Vector2f, Eigen::Vector2f>(vec1, vec2); }
  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::Vector2f)==vec2.type() ) { return InnerProduct<Eigen::VectorXf, Eigen::Vector2f>(vec1, vec2); }
  if( typeid(Eigen::Vector2f)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return InnerProduct<Eigen::Vector2f, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::Vector3f)==vec1.type() && typeid(Eigen::Vector3f)==vec2.type() ) { return InnerProduct<Eigen::Vector3f, Eigen::Vector3f>(vec1, vec2); }
  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::Vector3f)==vec2.type() ) { return InnerProduct<Eigen::VectorXf, Eigen::Vector3f>(vec1, vec2); }
  if( typeid(Eigen::Vector3f)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return InnerProduct<Eigen::Vector3f, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::Vector4f)==vec1.type() && typeid(Eigen::Vector4f)==vec2.type() ) { return InnerProduct<Eigen::Vector4f, Eigen::Vector4f>(vec1, vec2); }
  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::Vector4f)==vec2.type() ) { return InnerProduct<Eigen::VectorXf, Eigen::Vector4f>(vec1, vec2); }
  if( typeid(Eigen::Vector4f)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return InnerProduct<Eigen::Vector4f, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return InnerProduct<Eigen::VectorXf, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::Vector2i)==vec1.type() && typeid(Eigen::Vector2i)==vec2.type() ) { return InnerProduct<Eigen::Vector2i, Eigen::Vector2i>(vec1, vec2); }
  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::Vector2i)==vec2.type() ) { return InnerProduct<Eigen::VectorXi, Eigen::Vector2i>(vec1, vec2); }
  if( typeid(Eigen::Vector2i)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return InnerProduct<Eigen::Vector2i, Eigen::VectorXi>(vec1, vec2); }

  if( typeid(Eigen::Vector3i)==vec1.type() && typeid(Eigen::Vector3i)==vec2.type() ) { return InnerProduct<Eigen::Vector3i, Eigen::Vector3i>(vec1, vec2); }
  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::Vector3i)==vec2.type() ) { return InnerProduct<Eigen::VectorXi, Eigen::Vector3i>(vec1, vec2); }
  if( typeid(Eigen::Vector3i)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return InnerProduct<Eigen::Vector3i, Eigen::VectorXi>(vec1, vec2); }

  if( typeid(Eigen::Vector4i)==vec1.type() && typeid(Eigen::Vector4i)==vec2.type() ) { return InnerProduct<Eigen::Vector4i, Eigen::Vector4i>(vec1, vec2); }
  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::Vector4i)==vec2.type() ) { return InnerProduct<Eigen::VectorXi, Eigen::Vector4i>(vec1, vec2); }
  if( typeid(Eigen::Vector4i)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return InnerProduct<Eigen::Vector4i, Eigen::VectorXi>(vec1, vec2); }

  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return InnerProduct<Eigen::VectorXi, Eigen::VectorXi>(vec1, vec2); }

  // something went wrong
  assert(false);
  return 0.0;
}

boost::any EigenVectorAlgebra::OuterProduct(boost::any const& vec1, boost::any const& vec2) {
  if( typeid(Eigen::Vector2d)==vec1.type() && typeid(Eigen::Vector2d)==vec2.type() ) { return OuterProduct<Eigen::Matrix2d, Eigen::Vector2d, Eigen::Vector2d>(vec1, vec2); }
  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::Vector2d)==vec2.type() ) { return OuterProduct<Eigen::MatrixXd, Eigen::VectorXd, Eigen::Vector2d>(vec1, vec2); }
  if( typeid(Eigen::Vector2d)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return OuterProduct<Eigen::MatrixXd, Eigen::Vector2d, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::Vector3d)==vec1.type() && typeid(Eigen::Vector3d)==vec2.type() ) { return OuterProduct<Eigen::Matrix3d, Eigen::Vector3d, Eigen::Vector3d>(vec1, vec2); }
  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::Vector3d)==vec2.type() ) { return OuterProduct<Eigen::MatrixXd, Eigen::VectorXd, Eigen::Vector3d>(vec1, vec2); }
  if( typeid(Eigen::Vector3d)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return OuterProduct<Eigen::MatrixXd, Eigen::Vector3d, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::Vector4d)==vec1.type() && typeid(Eigen::Vector4d)==vec2.type() ) { return OuterProduct<Eigen::Matrix4d, Eigen::Vector4d, Eigen::Vector4d>(vec1, vec2); }
  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::Vector4d)==vec2.type() ) { return OuterProduct<Eigen::MatrixXd, Eigen::VectorXd, Eigen::Vector4d>(vec1, vec2); }
  if( typeid(Eigen::Vector4d)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return OuterProduct<Eigen::MatrixXd, Eigen::Vector4d, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::VectorXd)==vec1.type() && typeid(Eigen::VectorXd)==vec2.type() ) { return OuterProduct<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd>(vec1, vec2); }

  if( typeid(Eigen::Vector2f)==vec1.type() && typeid(Eigen::Vector2f)==vec2.type() ) { return OuterProduct<Eigen::Matrix2f, Eigen::Vector2f, Eigen::Vector2f>(vec1, vec2); }
  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::Vector2f)==vec2.type() ) { return OuterProduct<Eigen::MatrixXf, Eigen::VectorXf, Eigen::Vector2f>(vec1, vec2); }
  if( typeid(Eigen::Vector2f)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return OuterProduct<Eigen::MatrixXf, Eigen::Vector2f, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::Vector3f)==vec1.type() && typeid(Eigen::Vector3f)==vec2.type() ) { return OuterProduct<Eigen::Matrix3f, Eigen::Vector3f, Eigen::Vector3f>(vec1, vec2); }
  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::Vector3f)==vec2.type() ) { return OuterProduct<Eigen::MatrixXf, Eigen::VectorXf, Eigen::Vector3f>(vec1, vec2); }
  if( typeid(Eigen::Vector3f)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return OuterProduct<Eigen::MatrixXf, Eigen::Vector3f, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::Vector4f)==vec1.type() && typeid(Eigen::Vector4f)==vec2.type() ) { return OuterProduct<Eigen::Matrix4f, Eigen::Vector4f, Eigen::Vector4f>(vec1, vec2); }
  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::Vector4f)==vec2.type() ) { return OuterProduct<Eigen::MatrixXf, Eigen::VectorXf, Eigen::Vector4f>(vec1, vec2); }
  if( typeid(Eigen::Vector4f)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return OuterProduct<Eigen::MatrixXf, Eigen::Vector4f, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::VectorXf)==vec1.type() && typeid(Eigen::VectorXf)==vec2.type() ) { return OuterProduct<Eigen::MatrixXf, Eigen::VectorXf, Eigen::VectorXf>(vec1, vec2); }

  if( typeid(Eigen::Vector2i)==vec1.type() && typeid(Eigen::Vector2i)==vec2.type() ) { return OuterProduct<Eigen::Matrix2i, Eigen::Vector2i, Eigen::Vector2i>(vec1, vec2); }
  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::Vector2i)==vec2.type() ) { return OuterProduct<Eigen::MatrixXi, Eigen::VectorXi, Eigen::Vector2i>(vec1, vec2); }
  if( typeid(Eigen::Vector2i)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return OuterProduct<Eigen::MatrixXi, Eigen::Vector2i, Eigen::VectorXi>(vec1, vec2); }

  if( typeid(Eigen::Vector3i)==vec1.type() && typeid(Eigen::Vector3i)==vec2.type() ) { return OuterProduct<Eigen::Matrix3i, Eigen::Vector3i, Eigen::Vector3i>(vec1, vec2); }
  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::Vector3i)==vec2.type() ) { return OuterProduct<Eigen::MatrixXi, Eigen::VectorXi, Eigen::Vector3i>(vec1, vec2); }
  if( typeid(Eigen::Vector3i)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return OuterProduct<Eigen::MatrixXi, Eigen::Vector3i, Eigen::VectorXi>(vec1, vec2); }

  if( typeid(Eigen::Vector4i)==vec1.type() && typeid(Eigen::Vector4i)==vec2.type() ) { return OuterProduct<Eigen::Matrix4i, Eigen::Vector4i, Eigen::Vector4i>(vec1, vec2); }
  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::Vector4i)==vec2.type() ) { return OuterProduct<Eigen::MatrixXi, Eigen::VectorXi, Eigen::Vector4i>(vec1, vec2); }
  if( typeid(Eigen::Vector4i)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return OuterProduct<Eigen::MatrixXi, Eigen::Vector4i, Eigen::VectorXi>(vec1, vec2); }

  if( typeid(Eigen::VectorXi)==vec1.type() && typeid(Eigen::VectorXi)==vec2.type() ) { return OuterProduct<Eigen::MatrixXi, Eigen::VectorXi, Eigen::VectorXi>(vec1, vec2); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::AccessElement(boost::any const& vec, unsigned int const i) {
  if( typeid(Eigen::Vector2d)==vec.type() ) { return AccessElement<Eigen::Vector2d>(vec, i); }
  if( typeid(Eigen::Vector2f)==vec.type() ) { return AccessElement<Eigen::Vector2f>(vec, i); }
  if( typeid(Eigen::Vector2i)==vec.type() ) { return AccessElement<Eigen::Vector2i>(vec, i); }

  if( typeid(Eigen::Vector3d)==vec.type() ) { return AccessElement<Eigen::Vector3d>(vec, i); }
  if( typeid(Eigen::Vector3f)==vec.type() ) { return AccessElement<Eigen::Vector3f>(vec, i); }
  if( typeid(Eigen::Vector3i)==vec.type() ) { return AccessElement<Eigen::Vector3i>(vec, i); }

  if( typeid(Eigen::Vector4d)==vec.type() ) { return AccessElement<Eigen::Vector4d>(vec, i); }
  if( typeid(Eigen::Vector4f)==vec.type() ) { return AccessElement<Eigen::Vector4f>(vec, i); }
  if( typeid(Eigen::Vector4i)==vec.type() ) { return AccessElement<Eigen::Vector4i>(vec, i); }

  if( typeid(Eigen::VectorXd)==vec.type() ) { return AccessElement<Eigen::VectorXd>(vec, i); }
  if( typeid(Eigen::VectorXf)==vec.type() ) { return AccessElement<Eigen::VectorXf>(vec, i); }
  if( typeid(Eigen::VectorXi)==vec.type() ) { return AccessElement<Eigen::VectorXi>(vec, i); }

  // something went wront
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::Identity(std::type_info const& type, unsigned int const rows, unsigned int const cols) {
  if( type==typeid(Eigen::Vector2d) ) { return (Eigen::Matrix2d)Eigen::Matrix2d::Identity(); }
  if( type==typeid(Eigen::Vector2f) ) { return (Eigen::Matrix2f)Eigen::Matrix2f::Identity(); }
  if( type==typeid(Eigen::Vector2i) ) { return (Eigen::Matrix2i)Eigen::Matrix2i::Identity(); }
  
  if( type==typeid(Eigen::Vector3d) ) { return (Eigen::Matrix3d)Eigen::Matrix3d::Identity(); }
  if( type==typeid(Eigen::Vector3f) ) { return (Eigen::Matrix3f)Eigen::Matrix3f::Identity(); }
  if( type==typeid(Eigen::Vector3i) ) { return (Eigen::Matrix3i)Eigen::Matrix3i::Identity(); }
  
  if( type==typeid(Eigen::Vector4d) ) { return (Eigen::Matrix4d)Eigen::Matrix4d::Identity(); }
  if( type==typeid(Eigen::Vector4f) ) { return (Eigen::Matrix4f)Eigen::Matrix4f::Identity(); }
  if( type==typeid(Eigen::Vector4i) ) { return (Eigen::Matrix4i)Eigen::Matrix4i::Identity(); }

  if( type==typeid(Eigen::VectorXd) ) { return (Eigen::MatrixXd)Eigen::MatrixXd::Identity(rows, cols); }
  if( type==typeid(Eigen::VectorXf) ) { return (Eigen::MatrixXf)Eigen::MatrixXf::Identity(rows, cols); }
  if( type==typeid(Eigen::VectorXi) ) { return (Eigen::MatrixXi)Eigen::MatrixXi::Identity(rows, cols); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::Add(boost::any const& in0, boost::any const& in1) {
  // 2D vectors
  if( in0.type()==typeid(Eigen::Vector2d) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return Add<Eigen::Vector2d, Eigen::Vector2d>(in0, in1); }
    return Add<Eigen::Vector2d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2f) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return Add<Eigen::Vector2f, Eigen::Vector2f>(in0, in1); }
    return Add<Eigen::Vector2f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2i) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return Add<Eigen::Vector2i, Eigen::Vector2i>(in0, in1); }
    return Add<Eigen::Vector2i, Eigen::VectorXi>(in0, in1); 
  }

  // 3D vectors
  if( in0.type()==typeid(Eigen::Vector3d) ) {
    if( in1.type()==typeid(Eigen::Vector3d) ) { return Add<Eigen::Vector3d, Eigen::Vector3d>(in0, in1); }
    return Add<Eigen::Vector3d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3f) ) {
    if( in1.type()==typeid(Eigen::Vector3f) ) { return Add<Eigen::Vector3f, Eigen::Vector3f>(in0, in1); }
    return Add<Eigen::Vector3f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3i) ) {
    if( in1.type()==typeid(Eigen::Vector3i) ) { return Add<Eigen::Vector3i, Eigen::Vector3i>(in0, in1); }
    return Add<Eigen::Vector3i, Eigen::VectorXi>(in0, in1); 
  }

  // 4D vectors
  if( in0.type()==typeid(Eigen::Vector4d) ) {
    if( in1.type()==typeid(Eigen::Vector4d) ) { return Add<Eigen::Vector4d, Eigen::Vector4d>(in0, in1); }
    return Add<Eigen::Vector4d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4f) ) {
    if( in1.type()==typeid(Eigen::Vector4f) ) { return Add<Eigen::Vector4f, Eigen::Vector4f>(in0, in1); }
    return Add<Eigen::Vector4f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4i) ) {
    if( in1.type()==typeid(Eigen::Vector4i) ) { return Add<Eigen::Vector4i, Eigen::Vector4i>(in0, in1); }
    return Add<Eigen::Vector4i, Eigen::VectorXi>(in0, in1); 
  }

  // XD vectors
  if( in0.type()==typeid(Eigen::VectorXd) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return Add<Eigen::VectorXd, Eigen::Vector2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3d) ) { return Add<Eigen::VectorXd, Eigen::Vector3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4d) ) { return Add<Eigen::VectorXd, Eigen::Vector4d>(in0, in1); }
    return Add<Eigen::VectorXd, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXf) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return Add<Eigen::VectorXf, Eigen::Vector2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3f) ) { return Add<Eigen::VectorXf, Eigen::Vector3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4f) ) { return Add<Eigen::VectorXf, Eigen::Vector4f>(in0, in1); }
    return Add<Eigen::VectorXf, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXi) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return Add<Eigen::VectorXi, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return Add<Eigen::VectorXi, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return Add<Eigen::VectorXi, Eigen::Vector4i>(in0, in1); }
    return Add<Eigen::VectorXi, Eigen::VectorXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::Subtract(boost::any const& in0, boost::any const& in1) {
  // 2D vectors
  if( in0.type()==typeid(Eigen::Vector2d) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return Subtract<Eigen::Vector2d, Eigen::Vector2d>(in0, in1); }
    return Subtract<Eigen::Vector2d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2f) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return Subtract<Eigen::Vector2f, Eigen::Vector2f>(in0, in1); }
    return Subtract<Eigen::Vector2f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2i) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return Subtract<Eigen::Vector2i, Eigen::Vector2i>(in0, in1); }
    return Subtract<Eigen::Vector2i, Eigen::VectorXi>(in0, in1); 
  }

  // 3D vectors
  if( in0.type()==typeid(Eigen::Vector3d) ) {
    if( in1.type()==typeid(Eigen::Vector3d) ) { return Subtract<Eigen::Vector3d, Eigen::Vector3d>(in0, in1); }
    return Subtract<Eigen::Vector3d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3f) ) {
    if( in1.type()==typeid(Eigen::Vector3f) ) { return Subtract<Eigen::Vector3f, Eigen::Vector3f>(in0, in1); }
    return Subtract<Eigen::Vector3f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3i) ) {
    if( in1.type()==typeid(Eigen::Vector3i) ) { return Subtract<Eigen::Vector3i, Eigen::Vector3i>(in0, in1); }
    return Subtract<Eigen::Vector3i, Eigen::VectorXi>(in0, in1); 
  }

  // 4D vectors
  if( in0.type()==typeid(Eigen::Vector4d) ) {
    if( in1.type()==typeid(Eigen::Vector4d) ) { return Subtract<Eigen::Vector4d, Eigen::Vector4d>(in0, in1); }
    return Subtract<Eigen::Vector4d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4f) ) {
    if( in1.type()==typeid(Eigen::Vector4f) ) { return Subtract<Eigen::Vector4f, Eigen::Vector4f>(in0, in1); }
    return Subtract<Eigen::Vector4f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4i) ) {
    if( in1.type()==typeid(Eigen::Vector4i) ) { return Subtract<Eigen::Vector4i, Eigen::Vector4i>(in0, in1); }
    return Subtract<Eigen::Vector4i, Eigen::VectorXi>(in0, in1); 
  }

  // XD vectors
  if( in0.type()==typeid(Eigen::VectorXd) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return Subtract<Eigen::VectorXd, Eigen::Vector2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3d) ) { return Subtract<Eigen::VectorXd, Eigen::Vector3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4d) ) { return Subtract<Eigen::VectorXd, Eigen::Vector4d>(in0, in1); }
    return Subtract<Eigen::VectorXd, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXf) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return Subtract<Eigen::VectorXf, Eigen::Vector2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3f) ) { return Subtract<Eigen::VectorXf, Eigen::Vector3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4f) ) { return Subtract<Eigen::VectorXf, Eigen::Vector4f>(in0, in1); }
    return Subtract<Eigen::VectorXf, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXi) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return Subtract<Eigen::VectorXi, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return Subtract<Eigen::VectorXi, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return Subtract<Eigen::VectorXi, Eigen::Vector4i>(in0, in1); }
    return Subtract<Eigen::VectorXi, Eigen::VectorXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::ScalarMultiply(boost::any const& in0, boost::any const& in1) {
  if( in0.type()==typeid(double) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return ScalarMultiply<double, Eigen::Vector2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3d) ) { return ScalarMultiply<double, Eigen::Vector3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4d) ) { return ScalarMultiply<double, Eigen::Vector4d>(in0, in1); }
    
    return ScalarMultiply<double, Eigen::VectorXd>(in0, in1); 
  }

  if( in0.type()==typeid(float) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return ScalarMultiply<float, Eigen::Vector2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3f) ) { return ScalarMultiply<float, Eigen::Vector3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4f) ) { return ScalarMultiply<float, Eigen::Vector4f>(in0, in1); }

    return ScalarMultiply<float, Eigen::VectorXf>(in0, in1); 
  }

  if( in0.type()==typeid(int) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return ScalarMultiply<int, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return ScalarMultiply<int, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return ScalarMultiply<int, Eigen::Vector4i>(in0, in1); }

    return ScalarMultiply<int, Eigen::VectorXi>(in0, in1); 
  }

  if( in0.type()==typeid(unsigned int) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return ScalarMultiply<unsigned int, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return ScalarMultiply<unsigned int, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return ScalarMultiply<unsigned int, Eigen::Vector4i>(in0, in1); }

    return ScalarMultiply<unsigned int, Eigen::VectorXi>(in0, in1); 
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::Apply(boost::any const& A, boost::any const& x) {
  if( A.type()==typeid(Eigen::Vector2d) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return Apply<Eigen::Vector2d, Eigen::Vector2d>(A, x); }
    return Apply<Eigen::Vector2d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2f) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return Apply<Eigen::Vector2f, Eigen::Vector2f>(A, x); }
    return Apply<Eigen::Vector2f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2i) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return Apply<Eigen::Vector2i, Eigen::Vector2i>(A, x); }
    return Apply<Eigen::Vector2i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector3d) ) {
    if( x.type()==typeid(Eigen::Vector3d) ) { return Apply<Eigen::Vector3d, Eigen::Vector3d>(A, x); }
    return Apply<Eigen::Vector3d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3f) ) {
    if( x.type()==typeid(Eigen::Vector3f) ) { return Apply<Eigen::Vector3f, Eigen::Vector3f>(A, x); }
    return Apply<Eigen::Vector3f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3i) ) {
    if( x.type()==typeid(Eigen::Vector3i) ) { return Apply<Eigen::Vector3i, Eigen::Vector3i>(A, x); }
    return Apply<Eigen::Vector3i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector4d) ) {
    if( x.type()==typeid(Eigen::Vector4d) ) { return Apply<Eigen::Vector4d, Eigen::Vector4d>(A, x); }
    return Apply<Eigen::Vector4d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4f) ) {
    if( x.type()==typeid(Eigen::Vector4f) ) { return Apply<Eigen::Vector4f, Eigen::Vector4f>(A, x); }
    return Apply<Eigen::Vector4f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4i) ) {
    if( x.type()==typeid(Eigen::Vector4i) ) { return Apply<Eigen::Vector4i, Eigen::Vector4i>(A, x); }
    return Apply<Eigen::Vector4i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::VectorXd) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return Apply<Eigen::VectorXd, Eigen::Vector2d>(A, x); }
    if( x.type()==typeid(Eigen::Vector3d) ) { return Apply<Eigen::VectorXd, Eigen::Vector3d>(A, x); }
    if( x.type()==typeid(Eigen::Vector4d) ) { return Apply<Eigen::VectorXd, Eigen::Vector4d>(A, x); }

    return Apply<Eigen::VectorXd, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXf) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return Apply<Eigen::VectorXf, Eigen::Vector2f>(A, x); }
    if( x.type()==typeid(Eigen::Vector3f) ) { return Apply<Eigen::VectorXf, Eigen::Vector3f>(A, x); }
    if( x.type()==typeid(Eigen::Vector4f) ) { return Apply<Eigen::VectorXf, Eigen::Vector4f>(A, x); }

    return Apply<Eigen::VectorXf, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXi) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return Apply<Eigen::VectorXi, Eigen::Vector2i>(A, x); }
    if( x.type()==typeid(Eigen::Vector3i) ) { return Apply<Eigen::VectorXi, Eigen::Vector3i>(A, x); }
    if( x.type()==typeid(Eigen::Vector4i) ) { return Apply<Eigen::VectorXi, Eigen::Vector4i>(A, x); }

    return Apply<Eigen::VectorXi, Eigen::VectorXi>(A, x);
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::ApplyInverse(boost::any const& A, boost::any const& x) {
  if( A.type()==typeid(Eigen::Vector2d) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return ApplyInverse<Eigen::Vector2d, Eigen::Vector2d>(A, x); }
    return ApplyInverse<Eigen::Vector2d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2f) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return ApplyInverse<Eigen::Vector2f, Eigen::Vector2f>(A, x); }
    return ApplyInverse<Eigen::Vector2f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2i) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return ApplyInverse<Eigen::Vector2i, Eigen::Vector2i>(A, x); }
    return ApplyInverse<Eigen::Vector2i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector3d) ) {
    if( x.type()==typeid(Eigen::Vector3d) ) { return ApplyInverse<Eigen::Vector3d, Eigen::Vector3d>(A, x); }
    return ApplyInverse<Eigen::Vector3d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3f) ) {
    if( x.type()==typeid(Eigen::Vector3f) ) { return ApplyInverse<Eigen::Vector3f, Eigen::Vector3f>(A, x); }
    return ApplyInverse<Eigen::Vector3f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3i) ) {
    if( x.type()==typeid(Eigen::Vector3i) ) { return ApplyInverse<Eigen::Vector3i, Eigen::Vector3i>(A, x); }
    return ApplyInverse<Eigen::Vector3i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector4d) ) {
    if( x.type()==typeid(Eigen::Vector4d) ) { return ApplyInverse<Eigen::Vector4d, Eigen::Vector4d>(A, x); }
    return ApplyInverse<Eigen::Vector4d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4f) ) {
    if( x.type()==typeid(Eigen::Vector4f) ) { return ApplyInverse<Eigen::Vector4f, Eigen::Vector4f>(A, x); }
    return ApplyInverse<Eigen::Vector4f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4i) ) {
    if( x.type()==typeid(Eigen::Vector4i) ) { return ApplyInverse<Eigen::Vector4i, Eigen::Vector4i>(A, x); }
    return ApplyInverse<Eigen::Vector4i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::VectorXd) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return ApplyInverse<Eigen::VectorXd, Eigen::Vector2d>(A, x); }
    if( x.type()==typeid(Eigen::Vector3d) ) { return ApplyInverse<Eigen::VectorXd, Eigen::Vector3d>(A, x); }
    if( x.type()==typeid(Eigen::Vector4d) ) { return ApplyInverse<Eigen::VectorXd, Eigen::Vector4d>(A, x); }

    return ApplyInverse<Eigen::VectorXd, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXf) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return ApplyInverse<Eigen::VectorXf, Eigen::Vector2f>(A, x); }
    if( x.type()==typeid(Eigen::Vector3f) ) { return ApplyInverse<Eigen::VectorXf, Eigen::Vector3f>(A, x); }
    if( x.type()==typeid(Eigen::Vector4f) ) { return ApplyInverse<Eigen::VectorXf, Eigen::Vector4f>(A, x); }

    return ApplyInverse<Eigen::VectorXf, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXi) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return ApplyInverse<Eigen::VectorXi, Eigen::Vector2i>(A, x); }
    if( x.type()==typeid(Eigen::Vector3i) ) { return ApplyInverse<Eigen::VectorXi, Eigen::Vector3i>(A, x); }
    if( x.type()==typeid(Eigen::Vector4i) ) { return ApplyInverse<Eigen::VectorXi, Eigen::Vector4i>(A, x); }

    return ApplyInverse<Eigen::VectorXi, Eigen::VectorXi>(A, x);
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::Zero(std::type_info const& type, unsigned int const size) {  
  if( typeid(Eigen::Vector2d)==type ) { return (Eigen::Vector2d)Eigen::Vector2d::Zero(); }
  if( typeid(Eigen::Vector2f)==type ) { return (Eigen::Vector2f)Eigen::Vector2f::Zero(); }
  if( typeid(Eigen::Vector2i)==type ) { return (Eigen::Vector2i)Eigen::Vector2i::Zero(); }
  
  if( typeid(Eigen::Vector3d)==type ) { return (Eigen::Vector3d)Eigen::Vector3d::Zero(); }
  if( typeid(Eigen::Vector3f)==type ) { return (Eigen::Vector3f)Eigen::Vector3f::Zero(); }
  if( typeid(Eigen::Vector3i)==type ) { return (Eigen::Vector3i)Eigen::Vector3i::Zero(); }
  
  if( typeid(Eigen::Vector4d)==type ) { return (Eigen::Vector4d)Eigen::Vector4d::Zero(); }
  if( typeid(Eigen::Vector4f)==type ) { return (Eigen::Vector4f)Eigen::Vector4f::Zero(); }
  if( typeid(Eigen::Vector4i)==type ) { return (Eigen::Vector4i)Eigen::Vector4i::Zero(); }

  if( typeid(Eigen::VectorXd)==type ) { return (Eigen::VectorXd)Eigen::VectorXd::Zero(size); }
  if( typeid(Eigen::VectorXf)==type ) { return (Eigen::VectorXf)Eigen::VectorXf::Zero(size); }
  if( typeid(Eigen::VectorXi)==type ) { return (Eigen::VectorXi)Eigen::VectorXi::Zero(size); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any EigenVectorAlgebra::SquareRoot(boost::any const& obj) {
  if( typeid(Eigen::Vector2d)==obj.type() ) { return (Eigen::Vector2d)boost::any_cast<Eigen::Vector2d>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::Vector3d)==obj.type() ) { return (Eigen::Vector3d)boost::any_cast<Eigen::Vector3d>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::Vector4d)==obj.type() ) { return (Eigen::Vector4d)boost::any_cast<Eigen::Vector4d>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::VectorXd)==obj.type() ) { return (Eigen::VectorXd)boost::any_cast<Eigen::VectorXd>(obj).cwiseSqrt().matrix(); }

  if( typeid(Eigen::Vector2f)==obj.type() ) { return (Eigen::Vector2f)boost::any_cast<Eigen::Vector2f>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::Vector3f)==obj.type() ) { return (Eigen::Vector3f)boost::any_cast<Eigen::Vector3f>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::Vector4f)==obj.type() ) { return (Eigen::Vector4f)boost::any_cast<Eigen::Vector4f>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::VectorXf)==obj.type() ) { return (Eigen::VectorXf)boost::any_cast<Eigen::VectorXf>(obj).cwiseSqrt().matrix(); }

  if( typeid(Eigen::Vector2i)==obj.type() ) { return (Eigen::Vector2i)boost::any_cast<Eigen::Vector2i>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::Vector3i)==obj.type() ) { return (Eigen::Vector3i)boost::any_cast<Eigen::Vector3i>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::Vector4i)==obj.type() ) { return (Eigen::Vector4i)boost::any_cast<Eigen::Vector4i>(obj).cwiseSqrt().matrix(); }
  if( typeid(Eigen::VectorXi)==obj.type() ) { return (Eigen::VectorXi)boost::any_cast<Eigen::VectorXi>(obj).cwiseSqrt().matrix(); }

  // something went wrong
  assert(false);
  return boost::none;
}

double EigenVectorAlgebra::LogDeterminate(boost::any const& obj) {
  if( typeid(Eigen::Vector2d)==obj.type() ) { return LogDeterminate<Eigen::Vector2d>(obj); }
  if( typeid(Eigen::Vector2f)==obj.type() ) { return LogDeterminate<Eigen::Vector2f>(obj); }

  if( typeid(Eigen::Vector3d)==obj.type() ) { return LogDeterminate<Eigen::Vector3d>(obj); }
  if( typeid(Eigen::Vector3f)==obj.type() ) { return LogDeterminate<Eigen::Vector3f>(obj); }

  if( typeid(Eigen::Vector4d)==obj.type() ) { return LogDeterminate<Eigen::Vector4d>(obj); }
  if( typeid(Eigen::Vector4f)==obj.type() ) { return LogDeterminate<Eigen::Vector4f>(obj); }

  if( typeid(Eigen::VectorXd)==obj.type() ) { return LogDeterminate<Eigen::VectorXd>(obj); }
  if( typeid(Eigen::VectorXf)==obj.type() ) { return LogDeterminate<Eigen::VectorXf>(obj); }

  // something went wrong
  assert(false);
  return -1.0;
}
