#include "MUQ/Approximation/Polynomials/ProbabilistHermite.h"

#include <cmath>

using namespace muq::Approximation;


double ProbabilistHermite::DerivativeEvaluate(int const polyOrder, int const derivOrder, double const x) const {

    if((derivOrder > polyOrder) || (polyOrder==0))
        return 0.0;

    double c = 1.0;
    for(int k=polyOrder; k>polyOrder-derivOrder; --k)
        c *= k;

    return c*BasisEvaluate(polyOrder-derivOrder, x);

}

double ProbabilistHermite::Normalization(unsigned int polyOrder) const {
    return sqrt(2.0*M_PI) * std::tgamma(polyOrder+1);
}


double ProbabilistHermite::ak(unsigned int k) const {
  return 1.0;
}

double ProbabilistHermite::bk(unsigned int k) const {
  return 0.0;
}

double ProbabilistHermite::ck(unsigned int k) const {
  return k-1.0;
}


double ProbabilistHermite::phi0(double x) const {
  return 1.0;
}

double ProbabilistHermite::phi1(double x) const {
  return x;
}

REGISTER_SCALARBASIS_FAMILY(ProbabilistHermite)
