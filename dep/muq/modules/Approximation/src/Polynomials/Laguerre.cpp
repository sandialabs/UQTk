#include "MUQ/Approximation/Polynomials/Laguerre.h"

using namespace muq::Approximation;


double Laguerre::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = (derivOrder%2==0) ? 1.0 : -1.0;

    return c * Laguerre(a + derivOrder).BasisEvaluate(polyOrder-derivOrder,x);
}

double Laguerre::ak(unsigned int k) const{
  return -1.0/double(k);
}
double Laguerre::bk(unsigned int k) const{
  return (2.0*k+a-1.0)/double(k);
}
double Laguerre::ck(unsigned int k) const{
  return (k-1.0+a)/double(k);
}

double Laguerre::phi0(double x) const {
  return 1.0;
}

double Laguerre::phi1(double x) const {
  return 1.0 + a - x;
}

double Laguerre::Normalization(unsigned int polyOrder) const {
    return std::tgamma(polyOrder+a+1.0) / std::tgamma(polyOrder+1);
}


REGISTER_SCALARBASIS_FAMILY(Laguerre)
