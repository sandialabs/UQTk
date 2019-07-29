/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
//                          -*- C++ -*-

#ifndef _MyException_
#define _MyException_

#include <iostream>
#include <exception>
#include <string.h>

/**
 * Just an example exception - feel free to override this.
 */
class MyException : public std::exception {
public:
  /// Construct an exception using a C-style character string.
  MyException(const char* errormessage) {
    std::cerr << "ERROR:  " << errormessage << "\n";
    error_ = std::string("MyException: ") + errormessage;
  }

  /// Construct an exception using a C++-style string
  MyException(const std::string& errormessage) {
    std::cerr << "ERROR:  " << errormessage << "\n";
    error_ = std::string("MyException: ") + errormessage;
  }

  /// Destroy.
  virtual ~MyException() throw() {
  }

  /// What's going on?
  const char* what() const throw() {
    try {
      return error_.c_str();
    } catch(...) {
      ;/// This function is not permitted to throw exceptions.
    }
    return error_.c_str();
  }
  
private:
  std::string error_;
};

#endif // _MyException_
