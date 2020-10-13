#ifndef ANYALGEBRA2_H_
#define ANYALGEBRA2_H_

namespace muq {
  namespace Modeling {

      /**
    template<typename ObjectType>
    class AnyAlgebraImplementation
    {
    public:
        static unsigned int Size(ObjectType const& obj, int const dim=-1);

        static double Norm(ObjectType const& obj);
        static double InnerProduct(ObjectType const& obj1, ObjectType const& obj2);

        static ObjType OuterProduct(ObjectType const& obj1, ObjectType const& obj2);
        
    };
      */


      class AnyMat{

          AnyMat(std::shared_ptr<boost::any> objIn);
          
          unsigned Size(const int dim=-1);
          double Norm();

          AnyMat operator*(AnyMat const& otherMat) const;

          // Access operators
          double  operator(int i, int j=-1) const;
          double& operator(int i, int j=-1);

          AnyMat operator+(AnyMat const& otherMat) const;
          AnyMat operator-(AnyMat const& otherMat) const;

          AnyMat& operator+=(AnyMat const& otherMat);
          AnyMat& operator-=(AnyMat const& otherMat);

          AnyMat operator=(AnyMat const& otherMat);

          AnyMat Solve(AnyMat const& rhs);

          AnyMat Cholesky();
          
          AnyMat Zero(int rows, int cols=-1) const;
          AnyMat Ones(int rows, int cols=-1) const;

          bool IsZero() const;

          AnyMat Identity(int rows, int cols) const;

          double LogDeterminant() const;

          AnyMat Concatenate(AnyMat const& otherMat, int axis=-1) const;
          
      private:
          std::shared_ptr<boost::any> obj;
      };

      
    /// Implement a generic way to do algebric operations on boost::any's
    class AnyAlgebra2 {
    public:

      static unsigned int Size(boost::any const& obj, int const dim=-1);
      

      /// The norm of an object
      /**
	 @param[in] obj We need the norm of this object
	 \return The norm
       */
      double Norm(boost::any const& obj) const;

      /// The inner product between two vectors
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The inner product
       */
      double InnerProduct(boost::any const& vec1, boost::any const& vec2) const;

      /// The outer product between two vectors
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The outer product
       */
      boost::any OuterProduct(boost::any const& vec1, boost::any const& vec2) const;

      /// Access an element of a vector/matrix
      /**
	 The return type is whatever the elements of the vector/matrix are (doubles, ints, ect ...)
	 @param[in] obj The vector/matrix whose data we want to access
	 @param[in] i We want to access the \f$i^{th}\f$ element/row of the vector/matrix (defaults to 0)
	 @param[in] j We want to access the \f$j^{th}\f$ col of the matrix (defaults to 0)
	 \return The \f$(i,j)^{th}}\f$ element of the vector/matrix
       */
      boost::any AccessElement(boost::any const& obj, unsigned int const i = 0, unsigned int const j = 0) const;

      /// Compute a zero vector
      /** 
	  @param[in] type We need a zero object of this type
	  @param[in] rows The size of the vector (defaults to 0 because some types have implied sizes (e.g., Eigen::Vector2d)) or number of rows of the matrix
	  @param[in] cols The number of columns in the matrix (defaults to 0 but, again, some types imply a size)
       */
      boost::any Zero(std::type_info const& type, unsigned int const rows = 0, unsigned int const cols = 0) const;

      /// Determine if an object is the zero object
      /**
	 @param[in] obj An input object
	 \return true: if obj is the zero object, false: if obj is not the zero object
       */
      bool IsZero(boost::any const& obj) const;

      /// Compute an identity object
      /**
	 If the input type is a vector (e.g., Eigen::Vector), return an identity matrix of corresponding type (e.g., Eigen::Matrix)
	 @param[in] type The type---return an identity of this type
	 @param[in] rows The number of rows (e.g., for a matrix) defaults to 0 since some types imply the size
	 @param[in] cols The number of columns (e.g., for a matrix) defaults to 0 since some types imply the size
	 \return An identity of some type
       */
      boost::any Identity(std::type_info const& type, unsigned int const rows=0, unsigned int const cols=0) const;

      /// Compute an identity object
      /**
	 If the input type is a vector (e.g., Eigen::Vector), return an identity matrix of corresponding type (e.g., Eigen::Matrix)
	 @param[in] type The type---return an identity of this type
	 @param[in] rows The number of rows (e.g., for a matrix) defaults to 0 since some types imply the size
	 @param[in] cols The number of columns (e.g., for a matrix) defaults to 0 since some types imply the size
	 \return An identity of some type
       */
      boost::any Identity(std::type_index const& type, unsigned int const rows=0, unsigned int const cols=0) const;

      /// Add two objects together
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      boost::any Add(boost::any const& in0, boost::any const& in1) const;

      /// Subtract two objects 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      boost::any Subtract(boost::any const& in0, boost::any const& in1) const;

      /// Multiply two objects 
      /**
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0-in1)
       */
      boost::any Multiply(boost::any const& in0, boost::any const& in1) const;

      /// Apply the inverse of a matrix
      /**
	 If the input is a vector, treat is as the diagonal of a matrix
	 @param[in] A We are applying the inverse of this matrix
	 @param[in] x We are applying the inverse to this vector
	 \return The result \f$y=A^{-1} x\f$
       */
      boost::any ApplyInverse(boost::any const& A, boost::any const& x) const;

      /// Apply a matrix (mat-vec)
      /**
	 If the input is a vector, treat is as the diagonal of a matrix
	 @param[in] A We are applying this matrix
	 @param[in] x We are applying the matrix to this vector
	 \return The result \f$y=A x\f$
       */
      boost::any Apply(boost::any const& A, boost::any const& x) const;

      /// The inverse
      /**
	 @param[in] obj We need the inverse of this object
	 \return The inverse
       */
      boost::any Inverse(boost::any const& obj) const;

      /// Compute the square root of an object
      /**
	 In the vector case, compute the square root of each component.  In the matrix case, compute the Cholesky
	 @param[in] obj We need the square root of this object
	 \return The square root
       */
      boost::any SquareRoot(boost::any const& obj) const;

      /// Compute the log-determinate
      /**
	 In the vector case, compute the determinate of a diagonal matrix.
	 @param[in] obj We need the determinate of this object
	 \return The determinate
       */
      double LogDeterminate(boost::any const& obj) const;

      /// Combine to vectors into one
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The combined vector [vec1, vec2]
       */
      boost::any Concatenate(boost::any const& vec1, boost::any const& vec2) const;

    private:

      /// The size of an object (implemented by a child for non standard types)
      /**
	 For vectors/matrices, return the number of elements.
	 @param[in] obj We need the size of this object
	 \return The size
       */
      virtual unsigned int SizeImpl(boost::any const& obj) const;
      
      /// The norm of an object
      /**
	 @param[in] obj We need the norm of this object
	 \return The norm
       */
      virtual double NormImpl(boost::any const& obj) const;
            
      /// The inner product between two vectors
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The inner product
       */
      virtual double InnerProductImpl(boost::any const& vec1, boost::any const& vec2) const;

      /// The outer product between two vectors
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The outer product
       */
      virtual boost::any OuterProductImpl(boost::any const& vec1, boost::any const& vec2) const;

      /// Access an element of a vector
      /**
	 MUQ automatically checks for some common input types.  However, the user may need to overload this function for special types.
	 @param[in] vec The vector whose data we want to access
	 @param[in] i We want to access the \f$i^{th}\f$ element/row
	 \return The \f$i^{th}\f$ element/row of the vector
       */
      virtual boost::any AccessElementImpl(boost::any const& vec, unsigned int const i) const;

      /// Compute an identity object 
      /**
	 @param[in] type The type---return an identity of this type
	 @param[in] rows The number of rows (e.g., for a matrix)
	 @param[in] cols The number of columns (e.g., for a matrix) 
	 \return An identity of some type
       */
      virtual boost::any IdentityImpl(std::type_index const& type, unsigned int const rows, unsigned int const cols) const;
      
      /// Add two objects together
      /**
	 MUQ automatically checks for some common pairs.  However, the user may need to overload this function for special types.
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The addition of in0 and in1 (in0+in1)
       */
      virtual boost::any AddImpl(boost::any const& in0, boost::any const& in1) const;
      
      /// Subtract two objects 
      /**
	 MUQ automatically checks for some common pairs.  However, the user may need to overload this function for special types.
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The subtraction of in0 and in1 (in0-in1)
       */
      virtual boost::any SubtractImpl(boost::any const& in0, boost::any const& in1) const;

      /// Multiply two objects 
      /**
	 MUQ automatically checks for some common pairs.  However, the user may need to overload this function for special types.
	 @param[in] in0 The first input
	 @param[in] in1 The second input
	 \return The multiplication of in0 and in1 (in0*in1)
       */
      virtual boost::any MultiplyImpl(boost::any const& in0, boost::any const& in1) const;

      /// Apply a matrix (mat-vec)
      /**
	 If the input is a vector, treat is as the diagonal of a matrix
	 @param[in] A We are applying this matrix
	 @param[in] x We are applying the matrix to this vector
	 \return The result \f$y=A x\f$
       */
      virtual boost::any ApplyImpl(boost::any const& A, boost::any const& x) const;
      
      /// Apply the inverse of a matrix
      /**
	 If the input is a vector, treat is as the diagonal of a matrix
	 @param[in] A We are applying the inverse of this matrix
	 @param[in] x We are applying the inverse to this vector
	 \return The result \f$y=A^{-1} x\f$
       */
      virtual boost::any ApplyInverseImpl(boost::any const& A, boost::any const& x) const;

      /// Compute a zero object for boost::any
      /** 
	  @param[in] type We need a zero object of this type
	  @param[in] rows The size of the vector (defaults to 0 because some types have implied sizes (e.g., Eigen::Vector2d)) or number of rows of the matrix
	  @param[in] cols The number of columns in the matrix (defaults to 0 but, again, some types imply a size)
       */
      virtual boost::any ZeroImpl(std::type_info const& type, unsigned int const rows, unsigned int const cols) const;

      /// Determine if an object is the zero object
      /**
	 @param[in] obj An input object
	 \return true: if obj is the zero object, false: if obj is not the zero object
       */
      virtual bool IsZeroImpl(boost::any const& obj) const;

      /// The inverse
      /**
	 @param[in] obj We need the inverse of this object
	 \return The inverse
       */
      virtual boost::any InverseImpl(boost::any const& obj) const;

      /// Compute the square root of an object
      /**
	 In the vector case, compute the square root of each component.  In the matrix case, compute the Cholesky
	 @param[in] obj We need the square root of this object
	 \return The square root
       */
      virtual boost::any SquareRootImpl(boost::any const& obj) const;

      /// Compute the log-determinate
      /**
	 In the vector case, compute the determinate of a diagonal matrix.
	 @param[in] obj We need the determinate of this object
	 \return The determinate
       */
      virtual double LogDeterminateImpl(boost::any const& obj) const;

      /// Combine to vectors into one
      /**
	 @param[in] vec1 The first vector
	 @param[in] vec2 The second vector
	 \return The combined vector [vec1, vec2]
       */
      virtual boost::any ConcatenateImpl(boost::any const& vec1, boost::any const& vec2) const;

    };
  } // namespace Modeling
} // namespace muq

#endif
