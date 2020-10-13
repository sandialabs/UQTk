#include "MUQ/Modeling/LinearAlgebra/SundialsAlgebra.h"

using namespace muq::Modeling;

SundialsAlgebra::SundialsAlgebra() {}

SundialsAlgebra::~SundialsAlgebra() {}

#if MUQ_HAS_SUNDIALS==1
bool SundialsAlgebra::IsSundialsVector(std::type_info const& obj_type) {
  return typeid(N_Vector)==obj_type;
}
#endif

#if MUQ_HAS_SUNDIALS==1
unsigned int SundialsAlgebra::Size(boost::any const& vec) {
  const N_Vector& sun = boost::any_cast<N_Vector const&>(vec);

  return NV_LENGTH_S(sun);
}
#endif

#if MUQ_HAS_SUNDIALS==1
boost::any SundialsAlgebra::AccessElement(N_Vector const& vec, unsigned int const i) {
  // check the size
    assert(i<NV_LENGTH_S(vec));

    // return the ith element
    return NV_Ith_S(vec, i);
}
#endif
