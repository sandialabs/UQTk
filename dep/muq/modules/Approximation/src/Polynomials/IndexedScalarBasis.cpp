#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

#include "MUQ/Utilities/Exceptions.h"

#include <cstdlib>
#include <typeinfo>
#include <memory>
#include <cxxabi.h>

using namespace muq::Modeling;
using namespace muq::Approximation;

IndexedScalarBasis::IndexedScalarBasis() :
  WorkPiece(std::vector<std::string>({typeid(unsigned int).name(), typeid(double).name()}), // input types (order, point)
	    std::vector<std::string>({typeid(double).name()})) // output times (polynomial evaluation)
{}

void IndexedScalarBasis::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // extract the inputs
  const unsigned int order = boost::any_cast<unsigned int>(inputs[0]);
  const double x = boost::any_cast<double>(inputs[1]);

  // evaluate the polynomial
  outputs.resize(1);
  outputs[0] = BasisEvaluate(order, x);
}



std::shared_ptr<IndexedScalarBasis> IndexedScalarBasis::Construct(std::string const& polyName){

  auto map = GetScalarBasisMap();
  auto it = map->find(polyName);
  if(it == map->end()){
      std::string message = "The basis family, \"" + polyName + "\" has not been registered with the scalar basis factory.  Does the class exist?\n  Registered families include:\n";
      for(auto iter = map->begin(); iter!=map->end(); ++iter)
          message += "    " + iter->first + "\n";
      message += "\n";

      throw muq::NotRegisteredError(message);
      return nullptr;
  }else{
      return it->second();
  }
}

std::shared_ptr<IndexedScalarBasis::ScalarBasisMapType> IndexedScalarBasis::GetScalarBasisMap() {
  // define a static map from type to constructor
  static std::shared_ptr<ScalarBasisMapType> map;

  if( !map ) { // if the map has not yet been created ...
    // ... create the map
    map = std::make_shared<ScalarBasisMapType>();
  }

  return map;
}
