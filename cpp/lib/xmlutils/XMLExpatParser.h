/* =====================================================================================

                      The UQ Toolkit (UQTk) version 3.1.1
                          Copyright (2021) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright 2021 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

#ifndef _util_xml_class_XMLExpatParser_
#define _util_xml_class_XMLExpatParser_

#include "XMLParser.h"
#include <iostream>
#include <vector>
#include <expat.h>

/**
 * An XML parser that uses the Expat library to handle the gruntwork.
 * This class requires that the Expat be installed on your system.
 *
 * Expat is available <a href="http://www.jclark.com/xml/expat.html">
 * at the Expat site</a>
 *
 * This class may not be fully exception safe, since there is no
 * good way of enforcing that the Expat parser is destroyed cleanly.
 */
class XMLExpatParser : public XMLParser {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// Construct a new parser.
  XMLExpatParser();

private:
  /// Blocked copy constructor.  It is not safe to copy this object
  /// since the Expat parser may save state which cannot be duplicated.
  XMLExpatParser(const XMLExpatParser&);

public:
  /// Destructor.
  virtual ~XMLExpatParser() throw();

private:
  /// Blocked assignment operator.  Not for public consumption.
  XMLExpatParser& operator=(const XMLExpatParser&);

public:
  /// Parse the given input buffer and return a parse tree.
  RefPtr<XMLElement> parse(std::istream&);

private:
  /// The Expat parser.
  XML_Parser parser_;

  /// The path that we have traversed so far in building the tree.
  std::vector< RefPtr<XMLElement> > path_;

  /// The current leaf of the parse tree.
  RefPtr<XMLElement> leaf_;

  /// The method used to parse the start tag.
  void do_start(const XML_Char*, const XML_Char**);

  /// The method used to parse the end tag.
  void do_end(const XML_Char*);

  /// The method used to parse character (content) data.
  void do_character_data(const XML_Char*, int);

  /// Initialize the state of the parser.
  void init();

public:
  /// Static wrapper method used as a callback to get the 'start' tag.
  /// This method is for internal use only.  Calling this method
  /// directly will most likely result in a segmentation fault.
  static void start_(void*, const XML_Char*, const XML_Char**);

  /// Static wrapper method used as a callback to get the 'end' tag.
  /// This method is for internal use only.  Calling this method
  /// directly will most likely result in a segmentation fault.
  static void end_(void*, const XML_Char*);

  /// Static wrapper method used as a callback to get character data.
  /// This method is for internal use only.  Calling this method
  /// directly will most likely result in a segmentation fault.
  static void character_data_(void*, const XML_Char*, int);
};

#endif // _util_xml_class_XMLExpatParser_
