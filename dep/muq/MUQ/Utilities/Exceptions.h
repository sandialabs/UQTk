#ifndef MUQEXCEPTIONS_H
#define MUQEXCEPTIONS_H

#include <stdexcept>
#include <string>

namespace muq
{

    /** @defgroup Exceptions

     */
    /** @ingroup Exceptions
        @class NotImplementedError
        @brief Class for virtual base functions that are not implemented.
        @details In general, it's best to implement abstract class interfaces with pure virtual functions.  However,
                 there are some situations where not all children of the base class will implement a function.  This
                 exception is meant to be used in such a case.  It should be raised in the base classes virtual
                 function.  When children override this function, no exception will be thrown.
    */
    class NotImplementedError : public std::logic_error
    {
    public:
        NotImplementedError(std::string const& message) : std::logic_error(message){};

    };

    /** @ingroup Exceptions
        @class NotRegisteredError
        @brief Used when a child class has not been registered with the factory method.
        @details Throughout MUQ, we have many abstract parent classes that define general mathematical objects (like MCMC kernels or orthogonal polynomials).
                 We also use a class registration system to map strings to constructors for children of these abstract classes.  This error is thrown
                 when a string does not match any classes registered with a factory class.  For example, this exception will be thrown if "ClownShoes" is passed
                 to the Polynomial constructor.
    */
    class NotRegisteredError : public std::logic_error
    {
    public:
        NotRegisteredError(std::string const& message) : std::logic_error(message){};
    };

    /** @ingroup Exceptions
        @class ExternalLibraryError
        @brief A MUQ dependency has failed.
        @details Exception thrown when an external libraries returns a failed flag.
    */
    class ExternalLibraryError : public std::logic_error
    {
    public:
        ExternalLibraryError(std::string const& message) : std::logic_error(message){};
    };

    /** @class WrongSizeError
        @ingroup Exceptions
        @brief Exception to throw when matrices, vectors, or arrays are the wrong size.
    */
    class WrongSizeError : public std::length_error
    {
    public:
        WrongSizeError(std::string const& message) : std::length_error(message){};
    };


};



#endif // #ifndef MUQEXCEPTIONS
