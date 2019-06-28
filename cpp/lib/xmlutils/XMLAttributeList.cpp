/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0.4
                     Copyright (2017) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2017) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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
// The implementation of a -*-C++-*- attribute list for XML data.


#include "XMLAttributeList.h"
#include "MyException.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

// Utility routines to convert string to double.
inline double to_double(const std::string& value) {
  char* endptr;
  double retval = strtod(value.c_str(), &endptr);
  // Here you could do some error checking with endptr.
  // endptr should point one-past last parsed character.
  return retval;
}

// Utility routine to convert string to int.
inline int to_int(const std::string& value) {
  char* endptr;
  double retval = strtod(value.c_str(), &endptr);
  // Here you could do some error checking with endptr.
  // endptr should point one-past last parsed character.
  // You could also ensure that retval == int(retval)
  // or change the rounding.
  return int(retval);
}

// Utility routine to convert int to string.
// This would be much more classy using the stringstream classes,
// but those are (unfortunately) a little broken in some versions of gcc.
inline std::string to_string(int value) {
  char buffer[80];  // unlikely that we'll have more than 80 char int value.
  sprintf(buffer, "%i", value);
    return std::string(buffer);
}

// Utility routine to convert double to string.
// Same disclaimer applies as for to_string(int).
inline std::string to_string(double value) {
  char buffer[80];
  sprintf(buffer, "%f", value);
  return std::string(buffer);
}

//
// Default constructor.
//
XMLAttributeList::XMLAttributeList()
{}

//
// Blocked copy constructor.
//
XMLAttributeList::XMLAttributeList(const XMLAttributeList&) :
  Object()
{
  throw MyException("XMLAttributList: No copy constructor");
}

//
// Destructor.
//
XMLAttributeList::~XMLAttributeList() {
}

//
// Blocked assignment operator.
//
XMLAttributeList&
XMLAttributeList::operator=(const XMLAttributeList&) {
  throw MyException("XMLAttributeList: No assignment operator");
}

//
// Get the number of attributes in the list.
//
int XMLAttributeList::size() const {
  return attribute_.size();
}

//
// Return true if a key is defined.
//
bool XMLAttributeList::has(const std::string& key) const {
  return (get_location(key) != end());
}

//
// Return the value associated with the given key.
//
const std::string& XMLAttributeList::get(const std::string& key) const {
  const_iterator it = get_location(key);
  if(it == end())
    throw MyException
      (std::string("XMLAttributeList::get:  Invalid key \"") + key + "\"");
  else
    return it->second;
  }

//
// Return the value for the given key or the given default value.
//
std::string XMLAttributeList::get(const std::string& key,
				  const std::string& def) const
{
  const_iterator it = get_location(key);
  if(it == end())
    return def;
  else
    return it->second;
}

//
// Get the value for the given key as an integer value.
//
int XMLAttributeList::get_int(const std::string& key) const {
  const_iterator it = get_location(key);
  if(it == end()) {
    throw MyException
      (std::string("XMLAttributeList::get_int: Invalid key ") + key);
  }
  else {
    try {
      // Unfortunately, the lexical cast to integer is very fragile
      // (core dumps on some gcc versions when given a non-integer number).
      double temp = to_double(it->second);
      return int(temp);
    } catch(std::exception&) {
      throw MyException(std::string("XMLAttributeList::get_int: ") +
			it->second + " is not a valid integer");
    }
  }
}

//
// Get the value for the given key as an integer or return a defaul value.
//
int XMLAttributeList::get_int(const std::string& key, int def) const {
  const_iterator it = get_location(key);
  if(it == end()) {
    return def;
  }
  else {
    try {
      return to_int(it->second);
    } catch(std::exception&) {
      throw MyException(std::string("XMLAttributeList::get_int: \"") +
			it->second + "\" is not a valid integer");
    }
  }
}

//
// Get the value for the given key as a real value.
//
double XMLAttributeList::get_double(const std::string& key) const {
  const_iterator it = get_location(key);
  if(it == end()) {
    throw MyException
      (std::string("XMLAttributeList::get_double: Invalid key ") + key);
  }
  else {
    try {
      return to_double(it->second);
    } catch(std::exception&) {
      throw MyException(std::string("XMLAttributeList::get_double: \"")
			+ it->second + "\" is not a valid number");
    }
  }
}

//
// Get the value for the given key as an integer or return a defaul value.
//
double
XMLAttributeList::get_double(const std::string& key, double def) const
{
  const_iterator it = get_location(key);
  if(it == end()) {
    return def;
  }
  else {
    try {
      return to_double(it->second);
    } catch(std::exception&) {
      throw MyException(std::string("XMLAttributeList::get_double: \"")
			+ it->second + " is not a valid number");
    }
  }
}

//
// Get the value for the given key as a boolean value.
//
bool XMLAttributeList::get_bool(const std::string& key) const {
  const_iterator it = get_location(key);
  if(it == end())
    throw MyException
      (std::string("XMLAttributeList::get_bool: Invalid key ") + key);
  else
    return boolean_value(it->second, "get_bool");
}

//
// Get the value for the given key as a boolean value or return the default.
//
bool XMLAttributeList::get_bool(const std::string& key, bool def) const {
  const_iterator it = get_location(key);
  if(it == end())
    return def;
  else
    return boolean_value(it->second, "get_bool");
}

//
// Assign an attribute to the given key.
//
void XMLAttributeList::set(const std::string& key, const std::string& val)
{
  std::string low = key;
  //make_lower_case(low);
  attribute_[low] = val;
}

//
// Assign an integer value to the given key.
//
void XMLAttributeList::set_int(const std::string& key, int val) {
  std::string low = key;
  //make_lower_case(low);
  attribute_[low] = to_string(val);
}

//
// Assign a real value to the given key.
//
void XMLAttributeList::set_double(const std::string& key, double val) {
  std::string low = key;
  //make_lower_case(low);
  // lexical_cast seems to discard almost all my precision here
  // so I am using sprintf instead
  static const int capacity = 51;
  char buf[capacity + 1];
  buf[capacity] = '\0';
  // I would like to make this snprintf, but xlC doesn't like std::snprintf.
  sprintf(buf, "%.5e", val);
  attribute_[low] = buf;
}

//
// Assign a boolean value to the given key.
//
void XMLAttributeList::set_bool(const std::string& key, bool val) {
  std::string low = key;
  make_lower_case(low);
  if(val == true)
    attribute_[low] = "true";
  else
    attribute_[low] = "false";
}

//
// Get an iterator to the first element.
//
XMLAttributeList::iterator XMLAttributeList::begin() {
  return attribute_.begin();
}

//
// Get an iterator past the last element.
//
XMLAttributeList::iterator XMLAttributeList::end() {
  return attribute_.end();
}

//
// Get an iterator to the first element in a const context.
//
XMLAttributeList::const_iterator XMLAttributeList::begin() const{
  return attribute_.begin();
}

//
// Get an iterator past the last element in a const context.
//
XMLAttributeList::const_iterator XMLAttributeList::end() const {
  return attribute_.end();
}

//
// Private method to convert a std::string to lower case.
//
void XMLAttributeList::make_lower_case(std::string& str) const {
  const int length = str.size();
  for(int i = 0; i < length; ++i)
    str[i] = tolower(str[i]);
}

//
// Get the location of the given key.
// Returns end() if the key is not found.
//
XMLAttributeList::iterator
XMLAttributeList::get_location(const std::string& key)
{
  // make sure the key is lower case. (disabled by BD 3/13/08)
  std::string low = key;
  //make_lower_case(low);
  
  // find the string.
  return attribute_.find(low);
}

//
// Get the location of the given key in a const context.
// Returns end() if the key is not found.
//
XMLAttributeList::const_iterator
XMLAttributeList::get_location(const std::string& key) const
{
  // make sure the key is lower case. (disabled by BD 3/13/08)
  std::string low = key;
  //make_lower_case(low);
  
  // find the string.
  return attribute_.find(low);
}

//
// Return the boolean value of the given string.
//
bool XMLAttributeList::boolean_value(const std::string& str,
				     const char* where) const 
{
  // convert the string to lower case
  std::string low = str;
  make_lower_case(low);
  
  // test the string.
  if((low == "1") || (low == "true") || (low == "yes"))
    return true;
  else if((low == "0") || (low == "false") || (low == "no"))
    return false;
  else
    throw MyException(std::string("XMLAttributeList::") + where +
		      ": The string \"" + str + "\" is not a valid " +
		      "boolean value");
}

