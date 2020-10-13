#include "MUQ/Approximation/Polynomials/HermiteFunction.h"

#include <cmath>

using namespace muq::Approximation;

unsigned HermiteFunction::nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

double HermiteFunction::BasisEvaluate(int const order, double const x) const
{
  double scaling = std::pow( std::pow(2, order) *
                   std::tgamma(order+1) *
                   std::sqrt(M_PI), -0.5) * std::exp(-0.5*x*x);

  return scaling * polyBase->BasisEvaluate(order, x);
}


double HermiteFunction::DerivativeEvaluate(int const polyOrder,
                                           int const derivOrder,
                                           double const x) const
{
  if(derivOrder>polyOrder)
    return 0.0;

  if(derivOrder==1){
    return std::sqrt(0.5*polyOrder) * BasisEvaluate(polyOrder-1,x) +
           std::sqrt(0.5*(polyOrder+1)) * BasisEvaluate(polyOrder+1,x);
  }

  double result = 0.0;

  double nfact = std::tgamma(polyOrder + 1);

  for(int i=0; i<=derivOrder; ++i){
    result += nChoosek(derivOrder, i) * std::pow(2, 0.5*(derivOrder-i)) *
              std::sqrt(nfact / std::tgamma(polyOrder - derivOrder + i + 1)) *
              BasisEvaluate(polyOrder-derivOrder+i, x) * polyBase->BasisEvaluate(i, x);
  }

  return result;
}
