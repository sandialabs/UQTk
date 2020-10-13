#include "MUQ/Approximation/Polynomials/Legendre.h"

using namespace muq::Approximation;

Legendre::Legendre() : OrthogonalPolynomial() {}

Legendre::~Legendre() {}

double Legendre::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    if(derivOrder==1){
        return polyOrder / (x * x - 1.0) * (x * BasisEvaluate(polyOrder, x) - BasisEvaluate(polyOrder - 1, x));
    }else{
        // Use the fact that dp_{n+1}/dx = (2n+1) p_n + dp_{n-1} /dx
        return (2*(polyOrder-1) + 1) * DerivativeEvaluate(polyOrder-1, derivOrder-1, x) + DerivativeEvaluate(polyOrder-2, derivOrder, x);
    }
}

double Legendre::ak(unsigned int k) const {
  return (2.0*k-1.0) / k;
}

double Legendre::bk(unsigned int k) const {
  return 0.0;
}

double Legendre::ck(unsigned int k) const {
  return (k-1.0)/double(k);
}

double Legendre::phi0(double x) const {
  return 1.0;
}

double Legendre::phi1(double x) const {
  return x;
}

double Legendre::Normalization(unsigned int polyOrder) const {
    return 2.0/(2.0*polyOrder + 1.0);
}


REGISTER_SCALARBASIS_FAMILY(Legendre)
