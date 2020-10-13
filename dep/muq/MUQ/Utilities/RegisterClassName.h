#ifndef REGISTERCLASSNAME_H_
#define REGISTERCLASSNAME_H_

#include <memory>

//this file is based on a snippet: http://pastebin.com/dSTLt7vW

/**
  This file contains a basic prescription for performing string registration. It allows you to go from a string and get 
  a function pointer something that you can call like make_shared<T>.
  
  If you're trying to register a class, and were directed here, then add:
  REGISTER_YOURCLASSNAME(NAME) 
  where YOURCLASSNAME is the type of registration, and NAME is the unqualified name of the class to register.
  
  If you're creating a new type of registered class, follow the directions below.
*/

namespace muq {
  namespace Utilities {

      template<typename T>
      struct shared_factory{
	template<typename... Args>
        std::shared_ptr<T> operator()(Args... args){return std::make_shared<T>(args...);}
      };
    
  } // namespace Utilities
} // namespace muq

#endif
