
#ifndef EIGENTESTUTILS_H_
#define EIGENTESTUTILS_H_

#include <Eigen/Core>
#include <gtest/gtest.h>

//#include <MUQ/Utilities/EigenUtils.h>

///A gtest predicate for exact equality, used: EXPECT_PRED_FORMAT2(MatrixEqual,expected,actual)
namespace muq{
namespace Utilities{
        
    template<typename derivedA, typename derivedB>
    bool MatrixEqual(Eigen::MatrixBase<derivedA> const& a, Eigen::MatrixBase<derivedB> const& b)
    {
	//check rows and cols equal
	if ((a.cols() != b.cols()) || (a.rows() != b.rows())) {
	    return false;
	}
	
	//return true if they're all the same
	return (a.array() == b.array()).all();
    }
    
    template<typename derivedA, typename derivedB>
    bool MatrixApproxEqual(Eigen::MatrixBase<derivedA> const& a, Eigen::MatrixBase<derivedB> const& b, double const tol)
    {
	if ((a.cols() != b.cols()) || (a.rows() != b.rows())) { //check the length first
	    return false;
	}
	
	//return true if all the differences are less than the tol
	
	return ((a - b).array().abs() < tol).all();
    }
    
    template<typename derivedA, typename derivedB>
    ::testing::AssertionResult MatrixEqual(const char *m_expr, const char *n_expr, derivedA m, derivedB n)
    {
	if (MatrixEqual(m, n)) {
	    return ::testing::AssertionSuccess();
	}
	
	return ::testing::AssertionFailure() << "Value of: " << m_expr << " == " << n_expr << "\nExpected: \n" << m <<
	    "\nActual: \n" << n;
    }
    
///A gtest predicate for approximate equality, used: EXPECT_PRED_FORMAT3(MatrixApproxEqual,expected,actual,eps)
    template<typename derivedA, typename derivedB>
    ::testing::AssertionResult MatrixApproxEqual(const char *m_expr,
						 const char *n_expr,
						 const char *eps_expr,
						 derivedA    m,
						 derivedB    n,
						 double      eps)
    {
	if (MatrixApproxEqual(m, n, eps)) {
	    return ::testing::AssertionSuccess();
	}
	
	return ::testing::AssertionFailure() << "Value of: " << m_expr << " == " << n_expr << " to within eps = " << eps <<
	    "\nExpected: \n" << m << "\nActual: \n" << n;
    }

} // namespace Utilities
} // namespace muq

#endif //EIGENTESTUTILS_H_
