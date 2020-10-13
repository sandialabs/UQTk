#include "MUQ/Approximation/Polynomials/PhysicistHermite.h"

#include <cmath>

using namespace muq::Approximation;

PhysicistHermite::PhysicistHermite() : OrthogonalPolynomial() {}

PhysicistHermite::~PhysicistHermite() {}

double PhysicistHermite::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = 1.0;
    for(int k=polyOrder; k>polyOrder-derivOrder; --k)
        c *= 2.0*k;

    return c*BasisEvaluate(polyOrder-derivOrder, x);
}

double PhysicistHermite::ak(unsigned int k) const{
  return 2.0;
}
double PhysicistHermite::bk(unsigned int k) const{
  return 0.0;
}
double PhysicistHermite::ck(unsigned int k) const{
  return 2.0*(k-1.0);
}

double PhysicistHermite::Normalization(unsigned int polyOrder) const {
    return sqrt(M_PI) * pow(2.0, static_cast<double>(polyOrder)) * std::tgamma(polyOrder+1);
}

double PhysicistHermite::phi0(double x) const {
  return 1.0;
}

double PhysicistHermite::phi1(double x) const {
  return 2.0*x;
}

REGISTER_SCALARBASIS_FAMILY(PhysicistHermite)
