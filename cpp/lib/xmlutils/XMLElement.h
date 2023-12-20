/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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
//                                 -*- C++ -*-

#ifndef _util_xml_class_XMLElement_
#define _util_xml_class_XMLElement_

#include "Object.h"
#include "XMLAttributeList.h"
#include <vector>
#include <set>

/**
 * \class XMLElement
 *
 * This is the implementation of a node in an XML parse tree.
 * Each node contains the following three containers, any (or all) of which
 * may be empty.
 *
 * attributes (XMLAttributeList):  Contains attributes (key/value pairs).
 * children (Vector<XMLElement>):  Contains children of this node.
 * content (Vector<std::string>):  Contains content (text) data.
 * 
 * The implementation is very limited.  In particular, the following
 * advanced features are missing:
 *
 * encoding:  Support for character types other than char/std::string.
 * comments:  Allowing comment blocks to accompany each element.
 * other xml types (control statements, etc.).
 *
 * This implementation is heavily based on Kevin Long's XMLObject.
 */
class XMLElement : public Object {
  template <class T> friend class RefPtr;
  template <class T> friend class ConstRefPtr;
public:
  /// Construct a new xml element object and give it a label.
  XMLElement(const std::string&);
  
private:
  /// Blocked copy constructor.  
  /// \throw MyException.
  XMLElement(const XMLElement&);
  
public:
  /// Destructor.
  virtual ~XMLElement();
  
private:
  /// Blocked assignment operator.  
  /// \throw MyExcepiton.
  XMLElement& operator=(const XMLElement&);
  
public:
  /// Get the label of this node.
  const std::string& label() const;
  
  /// Assign a new label to this node.
  void set_label(const std::string&);
  
  /// Utility function to check how many attributes this element has.
  /// This amounts to the same as calling '.attributes().size()'
  int count_attributes() const;
  
  /// Get access to the attribute list.
  RefPtr<XMLAttributeList> attributes();
  
  /// Assign an attribute list to this element.
  void set_attributes(RefPtr<XMLAttributeList>);
  
  /// Utility function to check how many children this element has.
  int count_children() const;
  
  /// Get the child with the given index.
  /// \throw MyException if the index is invalid.
  RefPtr<XMLElement> get_child(int);
  
  /// Find the first instance of a child with a given
  /// label and return a pointer to it.
  /// \note Since child labels do not need to be unique, there may
  /// be multiple instances matching children
  /// \throw MyException if the child label can not be found
  /// \todo Make this more elegant with the STL find_if function
  RefPtr<XMLElement> get_child(const std::string&);
  
  /// Add a child to the back of the list.
  /// Ignored if the child is already in the list.
  /// \throw MyException if adding the child would
  /// result in a cyclic relationship.
  /// \throw MyException if the child holds a NULL pointer.
  void add_child(RefPtr<XMLElement>);
  
  /// Same as add_child, but this allows for repeating children
 void add_child_rpt(RefPtr<XMLElement>);
  
  /// Erase all child elements from this node.
  void clear_children();
  
  /// Utility function to check how many lines of text content
  /// are associated with this element.
  int count_content() const;
  
  /// Get a line of content by index.
  /// \throw MyException if the index is out of range.
  const std::string& get_content_line(int);
  
  /// Add a line of content.
  void add_content_line(const std::string&);
  
  /// Clear all text content.
  void clear_content();
  
  /// The iterator type returned for list of children
  //typedef std::vector< RefPtr<XMLElement> >::iterator child_iterator;
  
private:
  /// The label of this element.
  std::string label_;
  
  /// The list of attributes associated with this element.
  RefPtr<XMLAttributeList> attributes_;
  
  /// The list of children associated with this element.
  std::vector< RefPtr<XMLElement> > children_;
  
  /// The list of content associated with this element.
  std::vector<std::string> content_;
  
  /// A private routine called recursively to ensure that we don't
  /// have a cyclic relationship.
  void recurse(RefPtr<XMLElement>, std::set< RefPtr<XMLElement> >);
};

#endif // _util_xml_class_XMLElement_
