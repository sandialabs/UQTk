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
// The implementation part of a -*-C++-*- class
// that uses Expat to parse an XML file.
//
// Helgi
// March 11, 2002.

#include "XMLExpatParser.h"
#include "MyException.h"
#include <cstdio>
#include <iostream>

// Utility routine to convert int to string.
// This would be much more classy using the stringstream classes,
// but those are (unfortunately) a little broken in some versions of gcc.
inline std::string to_string(int value) {
  char buffer[80];  // unlikely that we'll have more than 80 char int value.
  sprintf(buffer, "%i", value);
    return std::string(buffer);
}

//
// Construct a new Expat parser.
//
XMLExpatParser::XMLExpatParser() :
  XMLParser(), parser_(NULL)
{
  init();
}

//
// Blocked copy constructor.
//
XMLExpatParser::XMLExpatParser(const XMLExpatParser&) :
  Object(), XMLParser()
{
  throw MyException("XMLExpatParser: No copy constructor.");
}

//
// Destructor.  Hopefully, this is exception safe enough.
//
XMLExpatParser::~XMLExpatParser() throw() {
  try {
    XML_ParserFree(parser_);    } catch(...) {
      ; // ignore exceptions at this point.
    }
}

//
// Blocked copy constructor.
//
XMLExpatParser& XMLExpatParser::operator=(const XMLExpatParser&) {
  throw MyException("XMLExpatParser: No assignment operator");
}

//
// Parse the given buffer and return a parse tree.
//
RefPtr<XMLElement> XMLExpatParser::parse(std::istream& buf) {
  //std::cerr << "XMLExpatParser::parse\n";
  
  int line_number = 0;
  bool done = false;
  std::string line;
  while(! done) {
    getline(buf, line);
    
    done = !(buf.good());
    ++line_number;
    
    //std::cerr << "XMLExpatParser:  line = \"" << line << "\ : done = " 
    //	<< done << "\n";
    
    if(XML_Parse(parser_, line.c_str(), line.length(), done) == 0) {
      // we only get here if we had a parse error
      std::cerr << "Parse error on \"" << line << "\":  done = "
		<< done << "\n";
      
      throw MyException
	(std::string("XMLExpatParser::parse:  \"") +
	 XML_ErrorString(XML_GetErrorCode(parser_)) +
	 "\" on line number " + to_string(line_number));
    }
  }
  
  // Return a parse tree and clear the data from this class.
  if(path_.size() < 1)
    throw MyException("XMLExpatParser::parse: No XML data found.");
  // I am not dealing with path_.size() > 1, which would indicate that
  // some tags were left unclosed.  These tags are automatically matched.
  
  RefPtr<XMLElement> ret = path_[0];
  path_.clear();
  init();
  
  return ret;
}

//
// Parse a start tag.
//
void XMLExpatParser::do_start(const XML_Char* lbl, const XML_Char** attr)
{
  //std::cerr << "XMLExpatParser::do_start \"" << lbl << "\"\n";
  
  // Construct a new element to contain this info.
  RefPtr<XMLElement> leaf = new XMLElement(lbl);
  for(int i = 0; attr[i] != NULL; i += 2)
    leaf->attributes()->set(attr[i], attr[i+1]);
  
  // Make sure we keep track of the path we have taken so far.
  if(path_.size() == 0) {
    path_.push_back(leaf);
  }
  else {
    // add this node as a child of the most recently parsed node.
    path_.back()->add_child(leaf);
    path_.push_back(leaf);
  }
}

//
// Parse an end tag.
//
void XMLExpatParser::do_end(const XML_Char* /*lbl*/) {
  //std::cerr << "XMLExpatParser::do_end \"" << lbl << "\"\n";
  
  // We don't really need to check that the tags match up (expat does that).
  if(path_.size() > 1) {
    // we don't want to pop the top node (which is when path_.size() == 1).
    path_.pop_back();
  }
}

//
// Parse character data (content).
// The content string is not NULL-terminated.
//
void XMLExpatParser::do_character_data(const XML_Char* data, int size) {
  //std::cerr << "XMLExpatParser::do_character_data\n";
  
  // This should never happen.
  if(path_.size() < 1 || size < 0) {
    throw MyException
      ("XMLExpatParser::parse:  Character data with no node.");
  }
  
  // Skip lines that contain only spaces and newline
  static const char* blank = " \t\n\r";
  std::string str(data, size);
  size_t first_char = str.find_first_not_of(blank);
  if(size != 0 && first_char != std::string::npos) {
    size_t last_char = str.find_last_not_of(blank);
    if(last_char != std::string::npos && last_char > first_char) {
      str = str.substr(first_char, last_char);
      path_.back()->add_content_line(str);
    }
  }
}

//
// Initialize the parser.
//
void XMLExpatParser::init() {
  // Unfortunately, the parser doesn't seem to clean up its state fully.
  if(parser_ != NULL)
    XML_ParserFree(parser_);
  
  // Build the parser
  parser_ = XML_ParserCreate(NULL);
  // Associate the parser with this object.
  XML_SetUserData(parser_, this);
  // Set up the 'start' and 'end' methods.
  XML_SetElementHandler(parser_,
			&XMLExpatParser::start_,
			&XMLExpatParser::end_);
  // Set up the character data handler.
  XML_SetCharacterDataHandler(parser_, &XMLExpatParser::character_data_);
}

//
// Static wrapper method to perform callback on 'start' tags.
//
void XMLExpatParser::start_(void* object,
			    const XML_Char* label,
			    const XML_Char** attributes)
{
  ( (XMLExpatParser*)object )->do_start(label, attributes);
}

//
// Static wrapper method to perform callback on 'end' tags.
//
void XMLExpatParser::end_(void* object, const XML_Char* label) {
  ( (XMLExpatParser*)object )->do_end(label);
}

//
// Static wrapper method to perform callback on character data.
//
void XMLExpatParser::character_data_(void* object,
				     const XML_Char* data,
				     int size)
{
  ( (XMLExpatParser*)object )->do_character_data(data, size);
}
