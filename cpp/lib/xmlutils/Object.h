/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.2
                          Copyright (2022) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
//                               -*- C++ -*-

#ifndef _base_class_Object_
#define _base_class_Object_

#include "RefPtr.h"

/**
 * \class Object
 * Base class for reference counted objects.
 *
 * Part of the Particle Simulation Toolkit (pst)
 *
 * The "friend" classes "RefPtr" and "ConstRefPtr" take care of the
 * reference counting and garbage collection.  This means that it 
 * should be safe to create an array of reference counted objects,
 * as long as you *do not* assign a reference counted pointer to 
 * any at the entries in the array at any time.
 */
class Object {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// Construct a new reference counted object with a zero reference count
  Object() : refs_(0) {
  }
  
  /// Destroy this object
  virtual ~Object() {
  }
  
  /// Returns the number of references that are held to this object
  long int reference_count() const {
    return refs_;
  }
  
protected:
  /// Enables the friends of the class to increment and decrement the
  /// reference count.
  long int reference_grab() const {
    return ++refs_;
  }
  
  long int reference_release() const {
    return --refs_;
  }    
  
private:
  mutable long int refs_;
};

#endif //_utility_ref_Object_
