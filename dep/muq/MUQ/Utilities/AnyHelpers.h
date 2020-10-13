#ifndef ANYHELPERS_H
#define ANYHELPERS_H

#include <boost/any.hpp>

#include "MUQ/config.h"

namespace muq
{
    namespace Utilities
    {
        /** @brief Class for easily casting boost::any's in assignment operations
@details  Consider the following code involving boost::any_cast
@code
double A = 1.0;
boost::any anyA = A;
double b = boost::any_cast<double>(anyA);
double const& aRef = boost::any_cast<double const&>(anyA);
@endcode
This class, AnyCast, attempts to streamline these casting operations and ensure
that references are used whenever possible.  Using AnyCast, the above code would
be replace by
@code
double A = 1.0;
boost::any anyA = A;
double b = AnyCast(anyA);
double& aRef = AnyCast(anyA);
@endcode
NOTE: This class stores a reference to the boost::any in question and errors may
      occur if an object of this class persists longer than the original boost::any.
      For this reason, we recommend using AnyCast in expressions like the code above;
        */
        class AnyCast
        {
        public:
          AnyCast(boost::any& objIn) : obj(objIn){};

          //template<typename T>
          //operator T(){ return boost::any_cast<T>(obj);};

          template<typename T>
          operator T&(){ return boost::any_cast<T&>(obj);};

        private:
          boost::any& obj;
        };

        /** The same as muq::Utilities::AnyCast, but using const references. */
        class AnyConstCast
        {
        public:
          AnyConstCast(boost::any const& objIn) : obj(objIn){};

          //template<typename T>
          //operator T(){ return boost::any_cast<T>(obj);};

          template<typename T>
          operator T const&(){ return boost::any_cast<T const&>(obj);};

#if MUQ_ANYCAST_COMPILES==0
          template<typename T>
	        operator T(){ return boost::any_cast<T const&>(obj);};
#endif

        private:
          boost::any const& obj;
        };
    } // namespace Utilities
} // namespace muq

#endif
