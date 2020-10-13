#ifndef SCALARALGEBRA_H_
#define SCALARALGEBRA_H_

#include <assert.h>
#include <cmath>

#include <boost/none.hpp>
#include <boost/any.hpp>

namespace muq {
  namespace Modeling {
    /// Linear algebra for scalar objects
    class ScalarAlgebra {
    public:

      ScalarAlgebra();

      virtual ~ScalarAlgebra();

      /// Is a boost::any a scalar type (double, float, int, or unsigned int)?
      /**
	 @param[in] obj_type We want to know if this object type is a scalar type
	 \return true: it is a scalar type, false: it is not a scalar type
       */
      static bool IsScalar(std::type_info const& obj_type);

      /// Determine if a scalar is zero 
      /**
	 @param[in] obj An input scalar
	 \return true: if obj is zero, false: if obj is not zero
       */
      static bool IsZero(boost::any const& obj);

      /// Compute a zero scalar
      /** 
	  @param[in] type We need a zero object of this type
       */
      static boost::any Zero(std::type_info const& type);

      /// Get the norm of a scalar (the magnitude)
      /**
	 @param[in] obj We need the magnitude of this scalar
	 \return The magnitude
       */
      static double Norm(boost::any const& obj);

      /// The inner product between two scalars
      /**
	 @param[in] vec1 The first scalar
	 @param[in] vec2 The second scalar
	 \return The inner product
       */
      static double InnerProduct(boost::any const& vec1, boost::any const& vec2);

      /// The outer product between two scalars
      /**
	 @param[in] vec1 The first scalar
	 @param[in] vec2 The second scalar
	 \return The outer product
       */
      static boost::any OuterProduct(boost::any const& vec1, boost::any const& vec2);

      /// Compute an identity object for a scalar
      /**
	 @param[in] type The type---return an identity of this type
	 \return An identity of some type
       */
      static boost::any Identity(std::type_info const& type);

      /// Add two scalars together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      static boost::any Add(boost::any const& in0, boost::any const& in1);

      /// Subtract two scalars
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      static boost::any Subtract(boost::any const& in0, boost::any const& in1);

      /// Multiply two scalars
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      static boost::any Multiply(boost::any const& in0, boost::any const& in1);

      /// The inverse
      /**
	 @param[in] obj We need the inverse of this object
	 \return The inverse
       */
      static boost::any Inverse(boost::any const& obj);

      /// Compute the square root of an object
      /**
	 @param[in] obj We need the square root of this object
	 \return The square root
       */
      static boost::any SquareRoot(boost::any const& obj);

      /// Compute the log-determinate
      /**
	 @param[in] obj We need the determinate of this object
	 \return The determinate
       */
      static double LogDeterminate(boost::any const& obj);

    private:

      /// The magnitude of a scalar
      /**
	 @param[in] obj We need the magnitude of this scalar
	 \return The magnitude
       */
      template<typename type>
	static inline double Magnitude(boost::any const& obj) {
	const double x = (double)boost::any_cast<type const>(obj);

	return std::abs(x);
      }

      /// Add two scalars together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      template<typename type0, typename type1> 
	static inline boost::any Add(boost::any const& in0, boost::any const& in1) {
	const type0 x0 = boost::any_cast<type0 const>(in0);
	const type1 x1 = boost::any_cast<type1 const>(in1);

	return x0 + x1;
      }

      /// Add two scalars together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      template<typename type0> 
	static inline boost::any Add(boost::any const& in0, boost::any const& in1) {
	if( in1.type()==typeid(double) ) { return Add<type0, double>(in0, in1); }
	if( in1.type()==typeid(float) ) { return Add<type0, float>(in0, in1); }
	if( in1.type()==typeid(int) ) { return Add<type0, int>(in0, in1); }
	if( in1.type()==typeid(unsigned int) ) { return Add<type0, unsigned int>(in0, in1); }

	// something went wrong
	assert(false);
	return boost::none;
      }

      /// Subtract two scalars
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      template<typename type0, typename type1> 
	static inline boost::any Subtract(boost::any const& in0, boost::any const& in1) {
	const type0 x0 = boost::any_cast<type0 const>(in0);
	const type1 x1 = boost::any_cast<type1 const>(in1);

	return x0 - x1;
      }

      /// Subtract two scalars 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      template<typename type0> 
	static inline boost::any Subtract(boost::any const& in0, boost::any const& in1) {
	if( in1.type()==typeid(double) ) { return Subtract<type0, double>(in0, in1); }
	if( in1.type()==typeid(float) ) { return Subtract<type0, float>(in0, in1); }
	if( in1.type()==typeid(int) ) { return Subtract<type0, int>(in0, in1); }
	if( in1.type()==typeid(unsigned int) ) { return Subtract<type0, unsigned int>(in0, in1); }

	// something went wrong
	assert(false);
	return boost::none;
      }

      /// Multiply two scalars
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      template<typename type0, typename type1> 
	static inline boost::any Multiply(boost::any const& in0, boost::any const& in1) {
	const type0 x0 = boost::any_cast<type0 const>(in0);
	const type1 x1 = boost::any_cast<type1 const>(in1);

	return x0 * x1;
      }

      /// Multiply two scalars 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      template<typename type0> 
	static inline boost::any Multiply(boost::any const& in0, boost::any const& in1) {
	if( in1.type()==typeid(double) ) { return Multiply<type0, double>(in0, in1); }
	if( in1.type()==typeid(float) ) { return Multiply<type0, float>(in0, in1); }
	if( in1.type()==typeid(int) ) { return Multiply<type0, int>(in0, in1); }
	if( in1.type()==typeid(unsigned int) ) { return Multiply<type0, unsigned int>(in0, in1); }

	// something went wrong
	assert(false);
	return boost::none;
      }
    };
  } // namespace Modeling
} // namespace muq

#endif
