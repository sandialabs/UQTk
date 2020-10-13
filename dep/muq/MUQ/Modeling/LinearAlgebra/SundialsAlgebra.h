#ifndef SUNDIALSALGEBRA_H_
#define SUNDIALSALGEBRA_H_

#include <assert.h>

#include <boost/none.hpp>
#include <boost/any.hpp>

#include "MUQ/config.h"

// Sundials includes
#if MUQ_HAS_SUNDIALS==1
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM
#endif

namespace muq {
  namespace Modeling {
    class SundialsAlgebra {
    public:
      SundialsAlgebra();

      ~SundialsAlgebra();

#if MUQ_HAS_SUNDIALS==1
      /// Is a boost::any an N_Vector type?
      /**
	 @param[in] obj_type We want to know if this object type is a N_Vector type
	 \return true: it is an N_Vector, false: it is not an N_Vector type
       */
      static bool IsSundialsVector(std::type_info const& obj);
#endif

#if MUQ_HAS_SUNDIALS==1
      /// The size of an N_Vector
      /**
	 @param[in] vec We will get the size of this vector
	 \return The size
       */
      static unsigned int Size(boost::any const& vec);
#endif
      
#if MUQ_HAS_SUNDIALS==1
      /// Access an element of a Sundials vector
      /**
	 @param[in] vec The vector whose data we want to access
	 @param[in] i We want to access the \f$i^{th}\f$ element of the vector
	 \return The \f$i^{th}\f$ element of the vector
       */
      static boost::any AccessElement(N_Vector const& obj, unsigned int const i);
#endif
      
    private:
    };
  }
} // namespace muq

#endif
