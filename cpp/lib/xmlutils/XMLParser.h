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
