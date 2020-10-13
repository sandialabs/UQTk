#include "MUQ/Approximation/Polynomials/Jacobi.h"

using namespace muq::Approximation;


double Jacobi::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = std::tgamma(a+b+polyOrder+1+derivOrder) / (std::pow(2.0, derivOrder) * std::tgamma(a+b+polyOrder+1));

    return c * Jacobi(a + derivOrder, b+derivOrder).BasisEvaluate(polyOrder-derivOrder, x);
}


double Jacobi::ak(unsigned int polyOrder) const{
  const double den = 2.0*polyOrder*(polyOrder + a + b);
  //return (2.0*polyOrder + a + b - 1.0)*(2.0*polyOrder + a + b) / den1;
  return (2.0*polyOrder+a+b-1)*(2.0*polyOrder+a+b)/den;
}
double Jacobi::bk(unsigned int polyOrder) const{
  const double den = 2.0*polyOrder*(polyOrder + a + b)*(2.0*polyOrder + a + b - 2.0);
  return (a*a-b*b) * (2.0*polyOrder + a + b - 1.0) / den;
}
double Jacobi::ck(unsigned int polyOrder) const{
  const double den = polyOrder*(polyOrder + a + b)*(2.0*polyOrder + a + b - 2.0);
  return (polyOrder + a-1.0)*(polyOrder + b-1.0)*(2.0*polyOrder + a + b) / den;
}


double Jacobi::phi0(double x) const {
    return 1.0;
}

double Jacobi::phi1(double x) const {
    return (0.5*(a+b)+1.0)*x + 0.5*(a-b);
}

double Jacobi::Normalization(unsigned int polyOrder) const {
    return (std::pow(2.0, a+b+1) / (2*polyOrder + a + b +1)) * std::tgamma(polyOrder+a+1)*std::tgamma(polyOrder+b+1)/(std::tgamma(polyOrder+a+b+1) * std::tgamma(polyOrder + 1));
}


REGISTER_SCALARBASIS_FAMILY(Jacobi)
