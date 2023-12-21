/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.4
                          Copyright (2023) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2023 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at https://github.com/sandialabs/UQTk/discussions
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
