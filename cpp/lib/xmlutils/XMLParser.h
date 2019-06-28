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
//                                 -*- C++ -*-

#ifndef _util_xml_class_XMLParser_
#define _util_xml_class_XMLParser_

#include "XMLElement.h"
#include <iostream>

/**
 * \class XMLParser
 * A pure abstract base class for parsers that read data from
 * an XML file and return the top node of a parse tree.
 * The parse tree node is a RefPtr< XMLElement >.
 */
class XMLParser : virtual public Object {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// Default constructor.  Intended for derived classes.
  XMLParser();
  
  /// Destructor.
  virtual ~XMLParser();
  
  /// Parse the given input buffer and return the parse tree.
  virtual RefPtr<XMLElement> parse(std::istream&) = 0;
};

#endif // _util_xml_class_XMLParser_
