#ifndef ANYWRITER_H
#define ANYWRITER_H

#include <boost/any.hpp>

namespace muq{
namespace Utilities{
    
    template<typename T>
    struct AnyWriter
    {
        template<typename DestType>
        void operator()(boost::any const& obj, DestType& dest){ dest = boost::any_cast<T>(obj); };
    };


} // namespace Utilities
} // namespace muq


#endif
