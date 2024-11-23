/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.5
                          Copyright (2024) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
//                               -*- C++ -*-

#ifndef _utility_ref_RefPtr_
#define _utility_ref_RefPtr_

#include "MyException.h"
#include <typeinfo>  // for dynamic_cast
#include <unistd.h>
#include <stddef.h>

/**
 * \class RefPtr
 * Reference counted pointer that gives the holder 
 * modification privileges to the pointee.
 *
 * Part of the Particle Simulation Toolkit (pst)
 */
template <class T>
class RefPtr {
public:
  /// Make the typename that this pointer holds accessible to other objects.
  typedef T Type;
  
  /// Construct a new RefPtr and initialize the pointee to NULL.
  RefPtr() : ptr_(NULL) {}
  
  /// Construct a new RefPtr and initialize the pointee as given.
  RefPtr(T* p) : ptr_(p) {
    grab();
  }
  
  /// Construct a new RefPtr and initialize to the given RefPtr pointee.
  RefPtr(const RefPtr<T>& p) : ptr_(p.ptr_) {
    grab();
  }
  
  /// Perform a static cast to initialize this pointee.
  /// This cast is only valid if T is a parent class of Other
  template <class Other>
  RefPtr(RefPtr<Other> p) : ptr_(static_cast<T*>(p.pointee())) {
    grab();
  }
  
  /// Destroy this RefPtr.
  ~RefPtr() {
    release();
  }
  
  /// Assign the value of this RefPtr to the given pointee.
  RefPtr<T>& operator=(T* p) {
    if(p != ptr_) {
      release();  // release the old pointer
      ptr_ = p;   // assign our value to this one
      grab();     // grab this pointer
    }
    return *this;
  }
  
  /// Assign the value of this RefPtr to the pointee of the given RefPtr.
  RefPtr<T>& operator=(const RefPtr<T>& p) {
    if(p.ptr_ != ptr_) {
      release();
      ptr_ = p.ptr_;
      grab();
    }
    return *this;
  }
  
  /// Use dynamic_cast to set the pointee to the
  /// pointer that was passed in, and return *this.
  /// The returned value is NULL if the cast fails.
  template <class Other>
  RefPtr<T>& cast(Other* p) {
    if(p != ptr_) {
      release();
      
      //std::cout << "DEBUG: Dynamic cast from type " << typeid(p).name()
      //	<< " to " << typeid(T).name() << std::endl; 
      
      ptr_ = dynamic_cast<T*>(p);
      if(p != NULL && ptr_ == NULL) {
	throw MyException
	  (std::string("RefPtr::cast(Other):  Failed dynamic cast from ")
	   + std::string(typeid(Other).name()) + std::string(" to ") +
	   std::string(typeid(Type).name()));
	
      }
      grab();
    }
    return *this;
  }
  
  /// Use dynamic_cast to set the pointee to the
  /// pointee of the RefPtr given, and return *this.
  /// The returned value is NULL if the cast fails.
  template <class Other>
  RefPtr<T>& cast(RefPtr<Other> p) {
    if(p.ptr_ != ptr_) {
      release();
      
      //std::cout << "DEBUG: Dynamic cast from type " 
      //	<< typeid(p.pointee()).name()
      //	<< " to " << typeid(T).name() << std::endl; 
      
      ptr_ = dynamic_cast<T*>(p.pointee());
      if(p != NULL && ptr_ == NULL) {
	throw MyException
	  (std::string("RefPtr::cast(Other):  Failed dynamic cast from ")
	   + std::string(typeid(Other).name()) + std::string(" to ") +
	   std::string(typeid(Type).name()));
      }
      grab();
    }
    return *this;
  }
  
  /// Return the pointee of this RefPtr.
  /// This will throw an exception if the pointee is NULL.
  T* operator->() const {
    if(ptr_ == NULL) {
      std::cerr << "RefPtr<" << typeid(T).name() 
		<< ">::operator->() const invoked on a null pointer\n";
      throw MyException("RefPtr::operator->() const");
    }
    return ptr_;
  }
  
  /// Return a reference to the pointee of this RefPtr.
  /// This will not work right if the pointee is NULL.
  T& operator*() const {
    if(ptr_ == NULL) {
      std::cerr << "RefPtr<" << typeid(T).name()
		<< ">::operator*() const invoked on a null pointer\n";
      throw MyException("RefPtr::operator*() const");
    } 
    return *ptr_;
  }
  
  /// Return the pointee of this RefPtr.
  T* pointee() {
    return ptr_;
  }
  
  /// Return the pointee of this RefPtr in a const context.
  const T* pointee() const {
    return ptr_;
  }
  
  /// Compare the pointee of this RefPtr with the given pointer.
  bool operator==(const T* p) const {
    return ptr_ == p;
  }
  
  /// Compare the value of this pointee with the pointee of the given RefPtr.
  bool operator==(const RefPtr<T>& p) const {
    return ptr_ == p.pointee();
  }
  
  /// Test inequality.
  bool operator!=(const T* p) const {
    return ptr_ != p;
  }
  
  /// Test inequality.
  bool operator!=(const RefPtr<T>& p) const {
    return ptr_ != p.ptr_;
  }
  
  /// Convenience routine to sort pointer values in standard containers.
  inline bool operator<(const RefPtr<T>& p) const {
    return ptr_ < p.ptr_;
  }
  
  /// Convenience routine to sort pointer values in standard containers.
  template <class Other>
  bool operator<(const RefPtr<Other>& p) const {
    return ptr_ < p.pointee();
  }
  
private:
  T* ptr_;
  
  /// Grab a reference to the current pointee if it is not NULL.
  inline void grab() {
    if(ptr_ != NULL) 
      ptr_->reference_grab();
  }
  
  /// Release the reference to the current pointee if it is not NULL.
  /// If this results in the reference count of the pointee dropping to zero,
  /// delete the object pointed to.
  inline void release() {
    if(ptr_ != NULL) {
      if(ptr_->reference_release() == 0)
	delete ptr_;
    }
  }
};

#endif //_utility_ref_RefPtr_
