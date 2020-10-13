#include "MUQ/Approximation/Quadrature/Quadrature.h"

#include "MUQ/Utilities/Exceptions.h"

using namespace muq::Approximation;

unsigned int Quadrature::Exactness(unsigned int quadOrder) const
{
  std::string className = typeid(*this).name();

  throw muq::NotImplementedError("The Exactness method has not been implmented for class \"" + className + "\".");
  return 0;
}
