#ifndef EIGENMATRIXALGEBRA_H_
#define EIGENMATRIXALGEBRA_H_

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <boost/none.hpp>
#include <boost/any.hpp>

namespace muq {
  namespace Modeling {
    /// Linear algebra for Eigen::Matrix's
    class EigenMatrixAlgebra {
    public:

      EigenMatrixAlgebra();

      ~EigenMatrixAlgebra();

      /// Is a boost::any an Eigen::Matrix type?
      /**
	 @param[in] obj_type We want to know if this object type is a Eigen::Matrix
	 \return true: it is an Eigen::Matrix2d, Eigen::Matrix3d, Eigen::Matrix4d, or Eigen::MatrixXd, false: it is not an Eigen::Matrixn type
       */
      static bool IsEigenMatrix(std::type_info const& obj_type);

      /// Determine if an Eigen::Matrix is zero 
      /**
	 @param[in] obj An input matrix
	 \return true: if obj is zero, false: if obj is not zero
       */
      static bool IsZero(boost::any const& obj);

      /// The norm of an Eigen::Matrix
      /**
	 @param[in] mat We will get the norm of this matrix
	 \return The norm
       */
      static double Norm(boost::any const& mat);

      /// The size of an Eigen::Matrix
      /**
	 @param[in] mat We will get the size of this matrix
	 @param[in] dim The dimension For matrices: for dim=-1 return the total number of elements, for dim=0, return the number of rows, for dim=1, return the number of colums
	 \return The size
       */
      static unsigned int Size(boost::any const& mat, int const dim);

      /// Access an element of an Eigen::Matrix
      /**
	 The return type is whatever the elements of the vector are (doubles, ints, ect ...)
	 @param[in] mat The matrix whose data we want to access
	 @param[in] i We want to access the \f$i^{th}\f$ row of the vector
	 @param[in] j We want to access the \f$j^{th}\f$ col of the vector
	 \return The \f$(i,j)^{th}\f$ element of the vector
       */
      static boost::any AccessElement(boost::any const& mat, unsigned int const i, unsigned int const j);

      /// Compute an identity Eigen::Matrix
      /**
	 @param[in] type The type---return an identity of this type
	 @param[in] rows The number of rows (e.g., for a matrix)
	 @param[in] cols The number of columns (e.g., for a matrix) 
	 \return An identity of some type
       */
      static boost::any Identity(std::type_info const& type, unsigned int const rows, unsigned int const cols);

      /// Add two Eigen::Matrices together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      static boost::any Add(boost::any const& in0, boost::any const& in1);

      /// Subtract two Eigen::Matrices
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0+in1)
       */
      static boost::any Subtract(boost::any const& in0, boost::any const& in1);

      /// Multiply two Eigen::Matrices together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0+in1)
       */
      static boost::any Multiply(boost::any const& in0, boost::any const& in1);

      /// Multiply Eigen::Matrix times scalars 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      static boost::any ScalarMultiply(boost::any const& in0, boost::any const& in1);

      /// Compute a zero Eigen::Matrix
      /** 
	  @param[in] type We need a zero object of this EigenMatrixType
	  @param[in] rows The number of rows
	  @param[in] cols The number of columns
       */
      static boost::any Zero(std::type_info const& type, unsigned int const rows, unsigned int const cols);

      /// Apply the inverse of a matrix 
      /**
	 @param[in] A We are applying the inverse of this matrix 
	 @param[in] x We are applying the inverse to this vector
	 \return The result \f$y=A^{-1} x\f$
       */
      static boost::any ApplyInverse(boost::any const& A, boost::any const& x);

      /// Apply a matrix 
      /**
	 @param[in] A We are applying this matrix 
	 @param[in] x We are applying the matrix to this vector
	 \return The result \f$y=A x\f$
       */
      static boost::any Apply(boost::any const& A, boost::any const& x);

      /// Compute the square root of an object
      /**
	 @param[in] obj We need the square root of this object
	 \return The square root
       */
      static boost::any SquareRoot(boost::any const& obj);

      /// Compute the log-determinate of a diagonal matrix
      /**
	 @param[in] obj We need the determinate of this object
	 \return The determinate
       */
      static double LogDeterminate(boost::any const& obj);

    private:

      /// Add two Eigen::Matrices together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      template<typename type0, typename type1>
	static inline boost::any Add(boost::any const& in0, boost::any const& in1) {
	const type0& x0 = boost::any_cast<type0 const&>(in0);
	const type1& x1 = boost::any_cast<type1 const&>(in1);
	assert(x0.rows()==x1.rows());
	assert(x0.cols()==x1.cols());

	return (type0)(x0+x1);
      }

      /// Add an Eigen::Matrice to an Eigen::Vector---treating the vector as a diagonal---
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      template<typename mattype, typename vectype>
	static inline boost::any AddVector(boost::any const& in0, boost::any const& in1) {
	mattype mat = boost::any_cast<mattype>(in0);
	const vectype& vec = boost::any_cast<vectype const&>(in1);
	assert(mat.rows()==vec.size());
	assert(mat.cols()==vec.size());
	
	//Eigen::VectorXd& diag = mat.diagonal();
	mat.diagonal() += vec;

	return mat;

	//Eigen:: result = mat+vec.asDiagonal();
	//return result;
      }

      /// The size of an Eigen::Matrix
      /**
	 @param[in] mat We will get the size of this matrix
	 @param[in] dim The dimension For matrices: for dim=-1 return the total number of elements, for dim=0, return the number of rows, for dim=1, return the number of colums
	 \return The size
       */
      template<typename type>
	static inline unsigned int Size(boost::any const& mat, int const dim) {
	const type& eig = boost::any_cast<type const&>(mat);
	switch( dim ) {
	case 0:
	  return eig.rows();
	case 1:
	  return eig.cols();
	default:
	  return eig.size();
	}
      }

      /// Determine if an Eigen::Matrix is zero 
      /**
	 @param[in] obj An input matrix
	 \return true: if obj is zero, false: if obj is not zero
       */
      template<typename EigenType>
	static inline bool IsZero(boost::any const& obj) {
	const EigenType& v = boost::any_cast<EigenType const&>(obj);

	return (v.array()==EigenType::Zero(v.rows(), v.cols()).array()).all();
      }

      /// Access an element of an Eigen::Matrix
      /**
	 The return type is whatever the elements of the vector are (doubles, ints, ect ...)
	 @param[in] mat The matrix whose data we want to access
	 @param[in] i We want to access the \f$i^{th}\f$ row of the vector
	 @param[in] j We want to access the \f$j^{th}\f$ col of the vector
	 \return The \f$(i,j)^{th}\f$ element of the vector
       */
      template<typename mattype>
	static inline boost::any AccessElement(boost::any const& mat, unsigned int const i, unsigned int const j) {
	// get a constant reference to the matrix
	const mattype& matref = boost::any_cast<const mattype&>(mat);
	
	// check the size
	assert(i<matref.rows());
	assert(j<matref.cols());
	
	// return ith element
	return matref(i,j);
      }
      
      /// Subtract two Eigen::Matrices
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      template<typename type0, typename type1>
	static inline boost::any Subtract(boost::any const& in0, boost::any const& in1) {
	const type0& x0 = boost::any_cast<type0 const&>(in0);
	const type1& x1 = boost::any_cast<type1 const&>(in1);
	assert(x0.rows()==x1.rows());
	assert(x0.cols()==x1.cols());

	return (type0)(x0-x1);
      }

      /// Multiply two Eigen::Matrices
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0-in1)
       */
      template<typename type0, typename type1>
	static inline boost::any Multiply(boost::any const& in0, boost::any const& in1) {
	const type0& x0 = boost::any_cast<type0 const&>(in0);
	const type1& x1 = boost::any_cast<type1 const&>(in1);
	assert(x0.cols()==x1.rows());

	return (type0)(x0*x1);
      }

      /// Multiply Eigen::Matrix times scalars 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      template<typename type0, typename type1> 
	static inline boost::any ScalarMultiply(boost::any const& in0, boost::any const& in1) {
	const type0& x0 = boost::any_cast<type0 const&>(in0);
	const type1& x1 = boost::any_cast<type1 const&>(in1);
	
	return (type1)(x0 * x1);
      }

      /// Apply the inverse of a matrix 
      /**
	 @param[in] A We are applying the inverse of this matrix 
	 @param[in] x We are applying the inverse to this vector
	 \return The result \f$y=A^{-1} x\f$
       */
      template<typename mattype, typename vectype>
	inline static boost::any ApplyInverse(boost::any const& A, boost::any const& x) {
	const mattype& mat = boost::any_cast<mattype const&>(A);
	const vectype& vec = boost::any_cast<vectype const&>(x);
	assert(mat.rows()==vec.size());
	assert(mat.cols()==vec.size());

	// solve the system
	const vectype soln = mat.colPivHouseholderQr().solve(vec);

	return soln;
      }

      /// Apply the inverse of a matrix 
      /**
	 @param[in] A We are applying the inverse of this matrix 
	 @param[in] x We are applying the inverse to this vector
	 \return The result \f$y=A^{-1} x\f$
       */
      template<typename mattype, typename vectype>
	inline static boost::any ApplyCholeskyInverse(boost::any const& A, boost::any const& x) {
	const Eigen::LLT<mattype>& chol = boost::any_cast<Eigen::LLT<mattype> const&>(A);
	const mattype& L = chol.matrixL();
  
	const vectype& vec = boost::any_cast<vectype const&>(x);
	assert(L.rows()==vec.size());
	assert(L.cols()==vec.size());
	
	// solve the system
	vectype soln = L.template triangularView<Eigen::Lower>().solve(vec);
	L.template triangularView<Eigen::Lower>().transpose().solveInPlace(soln);

	return soln;
      }

      /// Apply a matrix 
      /**
	 @param[in] A We are applying this matrix 
	 @param[in] x We are applying the matrix to this vector
	 \return The result \f$y=A x\f$
       */
      template<typename mattype, typename vectype>
	inline static boost::any Apply(boost::any const& A, boost::any const& x) {
	const mattype& mat = boost::any_cast<mattype const&>(A);
	const vectype& vec = boost::any_cast<vectype const&>(x);
	assert(mat.cols()==vec.size());

	return (vectype)(mat*vec);
      }

      /// Apply a matrix 
      /**
	 @param[in] A We are applying this matrix 
	 @param[in] x We are applying the matrix to this vector
	 \return The result \f$y=A x\f$
       */
      template<typename mattype, typename vectype>
	inline static boost::any ApplyCholesky(boost::any const& A, boost::any const& x) {
	const Eigen::LLT<mattype>& chol = boost::any_cast<Eigen::LLT<mattype> const&>(A);
	const mattype& mat = chol.matrixL();
	const vectype& vec = boost::any_cast<vectype const&>(x);
	assert(mat.cols()==vec.size());

	vectype soln = mat.template triangularView<Eigen::Lower>().transpose()*vec;
	soln = mat.template triangularView<Eigen::Lower>()*soln;

	return soln;
      }

      /// Compute the square root of an object
      /**
	 @param[in] obj We need the square root of this object
	 \return The square root
       */
      template<typename mattype>
	inline static boost::any SquareRoot(boost::any const& obj) {
	const Eigen::LLT<mattype>& mat = boost::any_cast<Eigen::LLT<mattype> const&>(obj);

	return (mattype)(mat.matrixL());
      }

      /// Compute the log-determinate of a diagonal matrix
      /**
	 @param[in] obj We need the determinate of this object
	 \return The determinate
       */
      template<typename mattype>
	inline static double LogDeterminate(boost::any const& obj) {
	const Eigen::LLT<mattype>& chol = boost::any_cast<Eigen::LLT<mattype> const&>(obj);
	const mattype& mat = chol.matrixL();

	return 2.0*mat.diagonal().array().log().sum();
      }
    };
  } // namespace Modeling
} // namespace muq

#endif
