#ifndef VARIADICMACROS_H
#define VARIADICMACROS_H


template<typename T> struct argument_type;
template<typename T, typename U> struct argument_type<T(U)> { typedef U type; };

/** In many cases, we have functions that take vectors as input, such as WorkPiece::EvaluateImpl, but we may not want to manually construct a vector very time.  For example,
its sometimes easier to call EvaluateImpl(a, b, c, d), instead of EvaluateImpl(std::vector<Type> list).   These macros define a series of functions using variadic templates
to concatenate a,b,c,d into a vector [a,b,c,d].
*/


#define STATIC_VARIADIC_TO_VECTOR_PART1(functionName, inputType, outputType) \
    template<typename... Args>                                         \
    inline static argument_type<void(outputType)>::type functionName(Args const&... args) {     \
        std::vector< argument_type<void(inputType)>::type> vec;                                    \
	return functionName(vec, args...);                             \
    }

#define STATIC_VARIADIC_TO_VECTOR_PART2(functionName, inputType, outputType) \
    template<typename... Args>                                                                                               \
    inline static  argument_type<void(outputType)>::type functionName(std::vector<argument_type<void(inputType)>::type>& bounds, argument_type<void(inputType)>::type const& ith, Args const&... args) {     \
        bounds.push_back(ith);                                                                                               \
	return functionName(bounds, args...);                                                                              \
    }                                                                                                                        \
    inline static  argument_type<void(outputType)>::type functionName(std::vector<argument_type<void(inputType)>::type>& bounds, argument_type<void(inputType)>::type const& last) {                         \
        bounds.push_back(last);                                                                                              \
        return functionName(bounds);                                                                                        \
    }

#define STATIC_VARIADIC_TO_VECTOR(functionName, inputType, outputType) \
    STATIC_VARIADIC_TO_VECTOR_PART1(functionName, inputType, outputType) \
    STATIC_VARIADIC_TO_VECTOR_PART2(functionName, inputType, outputType)

#define VARIADIC_TO_VECTOR(functionName, inputType, outputType) \
    template<typename... Args>                                         \
    inline argument_type<void(outputType)>::type functionName(Args const&... args) {     \
        std::vector< argument_type<void(inputType)>::type> vec;                                    \
	      return functionName(vec, args...);                             \
    }                                                                  \
    template<typename NextType, typename... Args>                                                                                               \
    inline argument_type<void(outputType)>::type functionName(std::vector<argument_type<void(inputType)>::type>& bounds, NextType const& ith, Args const&... args) {     \
        static_assert(std::is_same<argument_type<void(inputType)>::type, NextType>::value, "In "#functionName", cannot cast input to "#inputType"."); \
        bounds.push_back((argument_type<void(inputType)>::type&)ith);                                                                                               \
	      return functionName(bounds, args...);                                                                              \
    }                                                                                                                        \
    inline argument_type<void(outputType)>::type functionName(std::vector<argument_type<void(inputType)>::type>& bounds, argument_type<void(inputType)>::type const& last) {                         \
        bounds.push_back(last);                                                                                              \
        return functionName(bounds);                                                                                        \
    }

#define VARIADIC_TO_REFVECTOR(functionName, inputType, outputType) \
    template<typename... Args>                                         \
    inline argument_type<void(outputType)>::type functionName(Args const&... args) {     \
        ref_vector< argument_type<void(inputType)>::type> vec;                                    \
	      return functionName(vec, args...);                             \
    }                                                                  \
    template<typename NextType, typename... Args>                                                                                               \
    inline argument_type<void(outputType)>::type functionName(ref_vector<argument_type<void(inputType)>::type>& bounds, NextType const& ith, Args const&... args) {     \
        static_assert(std::is_same<argument_type<void(inputType)>::type, NextType>::value, "In "#functionName", cannot cast input to "#inputType"."); \
        bounds.push_back(std::cref((argument_type<void(inputType)>::type&)ith));                                                                                               \
	      return functionName(bounds, args...);                                                                              \
    }                                                                                                                        \
    inline argument_type<void(outputType)>::type functionName(ref_vector<argument_type<void(inputType)>::type>& bounds, argument_type<void(inputType)>::type const& last) {                         \
        bounds.push_back(std::cref(last));                                                                                              \
        return functionName(bounds);                                                                                        \
    }

#endif
