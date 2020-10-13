#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

#include "MUQ/Utilities/Exceptions.h"

#include <cstdlib>
#include <typeinfo>
#include <memory>
#include <cxxabi.h>

using namespace muq::Modeling;
using namespace muq::Approximation;


std::shared_ptr<OrthogonalPolynomial> OrthogonalPolynomial::Construct(std::string const& polyName)
{
  return std::dynamic_pointer_cast<OrthogonalPolynomial>(IndexedScalarBasis::Construct(polyName));
}


double OrthogonalPolynomial::Normalization(unsigned int polyOrder) const {

    std::string rawName = typeid(*this).name();

    int status = -4; // some arbitrary value to eliminate the compiler warning

    // enable c++11 by passing the flag -std=c++11 to g++
    std::unique_ptr<char, void(*)(void*)> res {
        abi::__cxa_demangle(rawName.c_str(), NULL, NULL, &status),
        std::free
    };

    std::string className = (status==0) ? res.get() : rawName;

    throw muq::NotImplementedError("The Normalization function has not been implemented for the class \"" + className + "\".  Is this polynomial family orthogonal?");

    return std::numeric_limits<double>::quiet_NaN();
};

double OrthogonalPolynomial::BasisEvaluate(int const order, double const x) const {

    if(order==0){
        return phi0(x);
    }else if(order==1){
        return phi1(x);
    }else{

      // "Downward" Clenshaw algorithm  http://mathworld.wolfram.com/ClenshawRecurrenceFormula.html
      double yk2 = 0.0;
      double yk1 = 0.0;
      double yk = 1.0;
      double alpha, beta;

      for( int k=order-1; k>=0; k-- ) {
        yk2 = yk1;
        yk1 = yk;

        alpha = ak(k+1)*x + bk(k+1);
        beta = -ck(k+2);
        yk = alpha*yk1 + beta*yk2;
      }
      beta = -ck(2);
      return yk1*phi1(x) + beta * phi0(x)*yk2;
    }
}

Eigen::VectorXd OrthogonalPolynomial::EvaluateAllTerms(int    const maxOrder,
                                                       double const x) const {

  Eigen::VectorXd output(maxOrder+1);
  output(0) = phi0(x);
  if(maxOrder>0)
    output(1) = phi1(x);

  for(int i=2; i<=maxOrder; ++i)
    output(i) = (ak(i)*x + bk(i))*output(i-1) - ck(i)*output(i-2);

  return output;
}
