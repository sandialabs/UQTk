#include "MUQ/Modeling/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Modeling;

AnyAlgebra::AnyAlgebra() {}

unsigned int AnyAlgebra::Size(boost::any const& obj, int const dim) const {
  // scalars are size one
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return 1; }

  // get the size of an Eigen::VectorXd
  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::Size(obj); }

  // get the size of an Eigen::MatrixXd
  if( EigenMatrixAlgebra::IsEigenMatrix(obj.type()) ) { return EigenMatrixAlgebra::Size(obj, dim); }

  // get the size of Sundials vectors
#if MUQ_HAS_SUNDIALS==1
  if( SundialsAlgebra::IsSundialsVector(obj.type()) ) { return SundialsAlgebra::Size(obj); }
#endif

  return SizeImpl(obj);
}

unsigned int AnyAlgebra::SizeImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: cannot compute the size of an object with type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::SizeImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::SizeImpl()" << std::endl << std::endl;

  return 0;
}

boost::any AnyAlgebra::Zero(std::type_info const& type, unsigned int rows, unsigned int const cols) const {
  if( ScalarAlgebra::IsScalar(type) ) { return ScalarAlgebra::Zero(type); }
  
  if( EigenVectorAlgebra::IsEigenVector(type) ) { return EigenVectorAlgebra::Zero(type, rows); }

  if( EigenMatrixAlgebra::IsEigenMatrix(type) ) { return EigenMatrixAlgebra::Zero(type, rows, cols); }

  if( type==typeid(double) ) {
    return 0.0;
  } 
  
  return ZeroImpl(type, rows, cols);
}

boost::any AnyAlgebra::ZeroImpl(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  std::cerr << std::endl << "ERROR: cannot compute zero of an object with type " << boost::core::demangle(type.name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ZeroImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ZeroImpl()" << std::endl << std::endl;

  return boost::none;
}

double AnyAlgebra::Norm(boost::any const& obj) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::Norm(obj); }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::Norm(obj); }

  if( EigenMatrixAlgebra::IsEigenMatrix(obj.type()) ) { return EigenMatrixAlgebra::Norm(obj); }
    
  return NormImpl(obj);
}

double AnyAlgebra::NormImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: Cannot compute the norm of an object with type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::NormImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::NormImpl()" << std::endl << std::endl;

  return -1.0;
}

double AnyAlgebra::InnerProduct(boost::any const& vec1, boost::any const& vec2) const {
  if( ScalarAlgebra::IsScalar(vec1.type()) && ScalarAlgebra::IsScalar(vec2.type()) ) { return ScalarAlgebra::InnerProduct(vec1, vec2); }

  if( EigenVectorAlgebra::IsEigenVector(vec1.type()) && EigenVectorAlgebra::IsEigenVector(vec2.type()) ) { return EigenVectorAlgebra::InnerProduct(vec1, vec2); }

  return InnerProductImpl(vec1, vec2);
}

double AnyAlgebra::InnerProductImpl(boost::any const& vec1, boost::any const& vec2) const {
  std::cerr << std::endl << "ERROR: Cannot compute the inner product between vectors with types " << boost::core::demangle(vec1.type().name()) << " and " << boost::core::demangle(vec2.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::InnerProductImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::InnerProductImpl()" << std::endl << std::endl;

  return 0.0;
}

boost::any AnyAlgebra::OuterProduct(boost::any const& vec1, boost::any const& vec2) const {
  if( ScalarAlgebra::IsScalar(vec1.type()) && ScalarAlgebra::IsScalar(vec2.type()) ) {
    return ScalarAlgebra::OuterProduct(vec1, vec2); }

  if( EigenVectorAlgebra::IsEigenVector(vec1.type()) && EigenVectorAlgebra::IsEigenVector(vec2.type()) ) { return EigenVectorAlgebra::OuterProduct(vec1, vec2); }

  return OuterProductImpl(vec1, vec2);
}

boost::any AnyAlgebra::OuterProductImpl(boost::any const& vec1, boost::any const& vec2) const {
  std::cerr << std::endl << "ERROR: Cannot compute the outer product between vectors with types " << boost::core::demangle(vec1.type().name()) << " and " << boost::core::demangle(vec2.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::OuterProductImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::OuterProductImpl()" << std::endl << std::endl;

  return boost::none;
}

bool AnyAlgebra::IsZero(boost::any const& obj) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::IsZero(obj); }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::IsZero(obj); }
  
  if( EigenMatrixAlgebra::IsEigenMatrix(obj.type()) ) { return EigenMatrixAlgebra::IsZero(obj); }
  
  return IsZeroImpl(obj);
}

bool AnyAlgebra::IsZeroImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: No way to determine if an object with type " << boost::core::demangle(obj.type().name()) << " is the zero vector." << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::IsZero()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::IsZero()" << std::endl << std::endl;

  return false;
}

boost::any AnyAlgebra::AccessElement(boost::any const& obj, unsigned int const i, unsigned int const j) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return obj; }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::AccessElement(obj, i); }

  if( EigenMatrixAlgebra::IsEigenMatrix(obj.type()) ) { return EigenMatrixAlgebra::AccessElement(obj, i, j); }

#if MUQ_HAS_SUNDIALS==1
  if( SundialsAlgebra::IsSundialsVector(obj.type()) ) { return SundialsAlgebra::AccessElement(boost::any_cast<const N_Vector&>(obj), i); }
#endif
  
  return AccessElementImpl(obj, i);
}

boost::any AnyAlgebra::AccessElementImpl(boost::any const& vec, unsigned int const i) const {
  std::cerr << std::endl << "ERROR: No way to access element " << i << " of a vector with type " << boost::core::demangle(vec.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::AccessElement()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::AccessElement()" << std::endl << std::endl;

  return boost::none;
}

boost::any AnyAlgebra::Identity(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  if( ScalarAlgebra::IsScalar(type) ) { return ScalarAlgebra::Identity(type); }

  if( EigenVectorAlgebra::IsEigenVector(type) ) { return EigenVectorAlgebra::Identity(type, rows, cols); }

  if( EigenMatrixAlgebra::IsEigenMatrix(type) ) { return EigenMatrixAlgebra::Identity(type, rows, cols); }
  
  return IdentityImpl(type, rows, cols);
}

boost::any AnyAlgebra::IdentityImpl(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  std::cerr << std::endl << "ERROR: No way to compute identy object with type " << boost::core::demangle(type.name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::IdentityImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::IdentityImpl()" << std::endl << std::endl;
  
  return boost::none;
}

boost::any AnyAlgebra::Add(boost::any const& in0, boost::any const& in1) const {
  // the first type is boost::none --- return the second
  if( in0.type()==typeid(boost::none) ) { return in1; }

  // the second type is boost::none --- return the first
  if( in1.type()==typeid(boost::none) ) { return in0; }
  
  if( ScalarAlgebra::IsScalar(in0.type()) && ScalarAlgebra::IsScalar(in1.type()) ) { return ScalarAlgebra::Add(in0, in1); }

  if( EigenVectorAlgebra::IsEigenVector(in0.type()) && EigenVectorAlgebra::IsEigenVector(in1.type()) ) { return EigenVectorAlgebra::Add(in0, in1); }

  if( EigenMatrixAlgebra::IsEigenMatrix(in0.type()) || EigenMatrixAlgebra::IsEigenMatrix(in1.type()) ) { return EigenMatrixAlgebra::Add(in0, in1); }
  
  return AddImpl(in0, in1);
}

boost::any AnyAlgebra::AddImpl(boost::any const& in0, boost::any const& in1) const {
  std::cerr << std::endl << "ERROR: No way to add type " << boost::core::demangle(in0.type().name()) << " and type " << boost::core::demangle(in1.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::AddImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::AddImpl()" << std::endl << std::endl;
  
  return boost::none;
}

boost::any AnyAlgebra::Subtract(boost::any const& in0, boost::any const& in1) const {
  if( ScalarAlgebra::IsScalar(in0.type()) && ScalarAlgebra::IsScalar(in1.type()) ) { return ScalarAlgebra::Subtract(in0, in1); }

  if( EigenVectorAlgebra::IsEigenVector(in0.type()) && EigenVectorAlgebra::IsEigenVector(in1.type()) ) { return EigenVectorAlgebra::Subtract(in0, in1); }

  if( EigenMatrixAlgebra::IsEigenMatrix(in0.type()) && EigenMatrixAlgebra::IsEigenMatrix(in1.type()) ) { return EigenMatrixAlgebra::Subtract(in0, in1); }
  
  // the first type is boost::none --- return the second
  if( in0.type()==typeid(boost::none) ) { return in1; }

  // the second type is boost::none --- return the first
  if( in1.type()==typeid(boost::none) ) { return in0; }

  return SubtractImpl(in0, in1);
}

boost::any AnyAlgebra::SubtractImpl(boost::any const& in0, boost::any const& in1) const {
  std::cerr << std::endl << "ERROR: No way to subtract type " << boost::core::demangle(in0.type().name()) << " and type " << boost::core::demangle(in1.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::SubtractImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::SubtractImpl()" << std::endl << std::endl;
  
  return boost::none;
}

boost::any AnyAlgebra::Multiply(boost::any const& in0, boost::any const& in1) const {
  if( ScalarAlgebra::IsScalar(in0.type()) && ScalarAlgebra::IsScalar(in1.type()) ) { return ScalarAlgebra::Multiply(in0, in1); }

  if( ScalarAlgebra::IsScalar(in0.type()) && EigenVectorAlgebra::IsEigenVector(in1.type()) ) { return EigenVectorAlgebra::ScalarMultiply(in0, in1); }
  if( ScalarAlgebra::IsScalar(in1.type()) && EigenVectorAlgebra::IsEigenVector(in0.type()) ) { return EigenVectorAlgebra::ScalarMultiply(in1, in0); }

  if( EigenMatrixAlgebra::IsEigenMatrix(in0.type()) && EigenMatrixAlgebra::IsEigenMatrix(in1.type()) ) { return EigenMatrixAlgebra::Multiply(in0, in1); }

  if( ScalarAlgebra::IsScalar(in0.type()) && EigenMatrixAlgebra::IsEigenMatrix(in1.type()) ) { return EigenMatrixAlgebra::ScalarMultiply(in0, in1); }
  if( ScalarAlgebra::IsScalar(in1.type()) && EigenMatrixAlgebra::IsEigenMatrix(in0.type()) ) { return EigenMatrixAlgebra::ScalarMultiply(in1, in0); }

  // the first type is boost::none --- return the second
  if( in0.type()==typeid(boost::none) ) { return in1; }

  // the second type is boost::none --- return the first
  if( in1.type()==typeid(boost::none) ) { return in0; }

  return MultiplyImpl(in0, in1);
}

boost::any AnyAlgebra::MultiplyImpl(boost::any const& in0, boost::any const& in1) const {
  std::cerr << std::endl << "ERROR: No way to multiply type " << boost::core::demangle(in0.type().name()) << " and type " << boost::core::demangle(in1.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::MultiplyImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::MultiplyImpl()" << std::endl << std::endl;
  
  return boost::none;
}

boost::any AnyAlgebra::ApplyInverse(boost::any const& A, boost::any const& x) const {
  if( ScalarAlgebra::IsScalar(A.type()) ) { return Multiply(Inverse(A), x); }

  if( EigenVectorAlgebra::IsEigenVector(A.type()) ) { return EigenVectorAlgebra::ApplyInverse(A, x); }
  
  if( EigenMatrixAlgebra::IsEigenMatrix(A.type()) ) { return EigenMatrixAlgebra::ApplyInverse(A, x); }
  
  return ApplyInverseImpl(A, x);
}

boost::any AnyAlgebra::ApplyInverseImpl(boost::any const& A, boost::any const& x) const {
  std::cerr << std::endl << "ERROR: No way to apply the inverse " << boost::core::demangle(A.type().name()) << " type to type " << boost::core::demangle(x.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ApplyInverseImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ApplyInverseImpl()" << std::endl << std::endl;
  
  return boost::none;
}

boost::any AnyAlgebra::Apply(boost::any const& A, boost::any const& x) const {
  if( ScalarAlgebra::IsScalar(A.type()) ) { return Multiply(A, x); }

  if( EigenVectorAlgebra::IsEigenVector(A.type()) ) { return EigenVectorAlgebra::Apply(A, x); }

  if( EigenMatrixAlgebra::IsEigenMatrix(A.type()) ) { return EigenMatrixAlgebra::Apply(A, x); }
  
  return ApplyImpl(A, x);
}

boost::any AnyAlgebra::ApplyImpl(boost::any const& A, boost::any const& x) const {
  std::cerr << std::endl << "ERROR: No way to apply " << boost::core::demangle(A.type().name()) << " type to type " << boost::core::demangle(x.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ApplyImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ApplyImpl()" << std::endl << std::endl;
  
  return boost::none;
}

boost::any AnyAlgebra::Inverse(boost::any const& obj) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::Inverse(obj); }
  
  return InverseImpl(obj);
}

boost::any AnyAlgebra::InverseImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: No way to compute the inverse of type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::InverseImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::InverseImpl()" << std::endl << std::endl;
  
  return boost::none;
}

boost::any AnyAlgebra::SquareRoot(boost::any const& obj) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::SquareRoot(obj); }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::SquareRoot(obj); }

  if( EigenMatrixAlgebra::IsEigenMatrix(obj.type()) ) { return EigenMatrixAlgebra::SquareRoot(obj); }
    
  return SquareRootImpl(obj);
}

boost::any AnyAlgebra::SquareRootImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: No way to compute the square root of type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::SquareRootImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::SquareRootImpl()" << std::endl << std::endl;
  
  return boost::none;
}

double AnyAlgebra::LogDeterminate(boost::any const& obj) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::LogDeterminate(obj); }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::LogDeterminate(obj); }

  if( EigenMatrixAlgebra::IsEigenMatrix(obj.type()) ) { return EigenMatrixAlgebra::LogDeterminate(obj); }
  
  return LogDeterminateImpl(obj);
}

double AnyAlgebra::LogDeterminateImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: No way to compute the determinate of type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::DeterminateImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::DeterminateImpl()" << std::endl << std::endl;
  
  return -1.0;
}

