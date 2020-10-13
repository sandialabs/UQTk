#ifndef EIGENVECTORALGEBRA_H_
#define EIGENVECTORALGEBRA_H_

#include <Eigen/Core>

#include <boost/none.hpp>
#include <boost/any.hpp>

namespace muq {
  namespace Modeling {
    /// Linear algebra for Eigen::Vector's
    class EigenVectorAlgebra {
    public:

      EigenVectorAlgebra();

      virtual ~EigenVectorAlgebra();

      /// Is a boost::any an Eigen::Vector type?
      /**
	 @param[in] obj_type We want to know if this object type is a Eigen::Vector type
	 \return true: it is an Eigen::Vector2d, Eigen::Vector3d, Eigen::Vector4d, or Eigen::VectorXd, false: it is not an Eigen::Vector type
       */
      static bool IsEigenVector(std::type_info const& obj_type);

      /// Determine if an Eigen::Vector is zero 
      /**
	 @param[in] obj An input vector
	 \return true: if obj is zero, false: if obj is not zero
       */
      static bool IsZero(boost::any const& obj);

      /// The size of an Eigen::Vector
      /**
	 @param[in] vec We will get the size of this vector
	 \return The size
       */
      static unsigned int Size(boost::any const& vec);

      /// The norm of an Eigen::Vector
      /**
	 @param[in] vec We will get the norm of this vector
	 \return The norm
       */
      static double Norm(boost::any const& vec);

      /// The inner product between two Eigen::Vector's
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The inner product
       */
      static double InnerProduct(boost::any const& vec1, boost::any const& vec2);

      /// The outer product between two Eigen::Vector's
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The outer product
       */
      static boost::any OuterProduct(boost::any const& vec1, boost::any const& vec2);

      /// Access an element of an Eigen::Vector
      /**
	 The return type is whatever the elements of the vector are (doubles, ints, ect ...)
	 @param[in] vec The vector whose data we want to access
	 @param[in] i We want to access the \f$i^{th}\f$ element of the vector
	 \return The \f$i^{th}\f$ element of the vector
       */
      static boost::any AccessElement(boost::any const& obj, unsigned int const i);

      /// Compute an identity Eigen::Matrix
      /**
	 @param[in] type The type---return an identity of this type
	 @param[in] rows The number of rows (e.g., for a matrix)
	 @param[in] cols The number of columns (e.g., for a matrix) 
	 \return An identity of some type
       */
      static boost::any Identity(std::type_info const& type, unsigned int const rows, unsigned int const cols);

      /// Add two Eigen::Vectors together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      static boost::any Add(boost::any const& in0, boost::any const& in1);

      /// Subtract two Eigen::Vectors 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      static boost::any Subtract(boost::any const& in0, boost::any const& in1);

      /// Multiply Eigen::Vector times scalars 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      static boost::any ScalarMultiply(boost::any const& in0, boost::any const& in1);

      /// Apply a diagonal matrix---the input vector is the diaganal of the matrix
      /**
	 @param[in] A We are applying this matrix
	 @param[in] x We are applying the matrix to this vector
	 \return The result \f$y=A x\f$
       */
      static boost::any Apply(boost::any const& A, boost::any const& x);

      /// Apply the inverse of a diagonal matrix 
      /**
	 @param[in] A We are applying the inverse of a matrix with this diagonal
	 @param[in] x We are applying the inverse to this vector
	 \return The result \f$y=A^{-1} x\f$
       */
      static boost::any ApplyInverse(boost::any const& A, boost::any const& x);

      /// Compute a zero Eigen::Vector
      /** 
	  @param[in] type We need a zero object of this EigenVectorType
	  @param[in] size The size of the vector 
       */
      static boost::any Zero(std::type_info const& type, unsigned int const size);

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
      
      /// Determine if an Eigen::Vector is zero 
      /**
	 @param[in] obj An input vector
	 \return true: if obj is zero, false: if obj is not zero
       */
      template<typename EigenType>
	static inline bool IsZero(boost::any const& obj) {
	const EigenType& v = boost::any_cast<EigenType const&>(obj);

	return (v.array()==EigenType::Zero(v.size()).array()).all();
      }
      
      /// The inner product between two Eigen::Vector's
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The inner product
       */
      template<typename EigenType1, typename EigenType2>
	static inline double InnerProduct(boost::any const& vec1, boost::any const& vec2) {
	const EigenType1& x1 = boost::any_cast<EigenType1 const&>(vec1);
	const EigenType2& x2 = boost::any_cast<EigenType2 const&>(vec2);
	assert(x1.size()==x2.size());
	
	return x1.dot(x2);    
      }

      /// The inner product between two Eigen::Vector's
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The inner product
       */
      template<typename mattype, typename EigenType1, typename EigenType2>
	static inline boost::any OuterProduct(boost::any const& vec1, boost::any const& vec2) {
	const EigenType1& x1 = boost::any_cast<EigenType1 const&>(vec1);
	const EigenType2& x2 = boost::any_cast<EigenType2 const&>(vec2);
	
	return (mattype)(x1*x2.transpose());    
      }

      /// Access an element of an Eigen::Vector
      /**
	 The return type is whatever the elements of the vector are (doubles, ints, ect ...)
	 @param[in] vec The vector whose data we want to access
	 @param[in] i We want to access the \f$i^{th}\f$ element of the vector
	 \return The \f$i^{th}\f$ element of the vector
       */
      template<typename vectype>
	static inline boost::any AccessElement(boost::any const& vec, unsigned int const i) {
	// get a constant reference to the vector
	const vectype& vecref = boost::any_cast<const vectype&>(vec);
	
	// check the size
	assert(i<vecref.size());
	
	// return ith element
	return vecref(i);
      }

      /// Add two Eigen::Vectors together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      template<typename type0, typename type1>
	static inline boost::any Add(boost::any const& in0, boost::any const& in1) {
	const type0& x0 = boost::any_cast<type0 const&>(in0);
	const type1& x1 = boost::any_cast<type1 const&>(in1);
	assert(x0.size()==x1.size());

	return (type0)(x0+x1);
      }

      /// Subtract two Eigen::Vectors
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      template<typename type0, typename type1>
	static inline boost::any Subtract(boost::any const& in0, boost::any const& in1) {
	const type0& x0 = boost::any_cast<type0 const&>(in0);
	const type1& x1 = boost::any_cast<type1 const&>(in1);
	assert(x0.size()==x1.size());

	return (type0)(x0-x1);
      }

      /// Multiply Eigen::Matrix times scalars 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      template<typename type0, typename type1> 
	static inline boost::any ScalarMultiply(boost::any const& in0, boost::any const& in1) {
	const type0 x0 = boost::any_cast<type0 const>(in0);
	const type1 x1 = boost::any_cast<type1 const>(in1);
	
	return (type1)(x0 * x1);
      }

      /// Apply a diagonal matrix---the input vector is the diaganal of the matrix
      /**
	 @param[in] A We are applying this matrix
	 @param[in] x We are applying the matrix to this vector
	 \return The result \f$y=A x\f$
       */
      template<typename mattype, typename vectype>
	static inline boost::any Apply(boost::any const& A, boost::any const& x) {
	const mattype& mat = boost::any_cast<mattype const&>(A);
	const vectype& vec = boost::any_cast<vectype const&>(x);
	assert(mat.size()==vec.size());
	
	return (vectype)(mat.asDiagonal()*vec);
      }

      /// Apply the inverse of a diagonal matrix 
      /**
	 @param[in] A We are applying the inverse of a matrix with this diagonal
	 @param[in] x We are applying the inverse to this vector
	 \return The result \f$y=A^{-1} x\f$
       */
      template<typename mattype, typename vectype>
	static inline boost::any ApplyInverse(boost::any const& A, boost::any const& x) {
	const mattype& mat = boost::any_cast<mattype const&>(A);
	const vectype& vec = boost::any_cast<vectype const&>(x);
	assert(mat.size()==vec.size());
	
	return (vectype)((1/mat.array()).matrix().asDiagonal()*vec);
      }
      
      /// Compute the log-determinate of a diagonal matrix
      /**
	 @param[in] obj We need the determinate of this object
	 \return The determinate
       */
      template<typename type>
	inline static double LogDeterminate(boost::any const& obj) {
	const type& mat = boost::any_cast<type const&>(obj);

	return mat.array().log().sum();
      }
    };
  } // namespace Modeling
} // namespace muq

#endif
