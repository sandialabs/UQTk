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
// The implementation of a -*-C++-*- class to hold xml element data.
//
// Helgi
// Feb. 25, 2002.

#include "XMLElement.h"
#include "MyException.h"
#include <algorithm>  // std::find

  //
  // Construct a new node and give it a label.
  //
  XMLElement::XMLElement(const std::string& lbl) :
    Object(), label_(lbl), attributes_(new XMLAttributeList)
  {}

  //
  // Blocked copy constructor.
  //
  XMLElement::XMLElement(const XMLElement&) :
    Object()
  {
    throw MyException("XMLElement: Blocked copy constructor");
  }

  //
  // Destructor.
  //
  XMLElement::~XMLElement() {
  }

  //
  // Blocked assignment operator.
  //
  XMLElement& XMLElement::operator=(const XMLElement&) {
    throw MyException("XMLElement: Blocked assignment operator");
  }

  //
  // Get the label of this node.
  //
  const std::string& XMLElement::label() const {
    return label_;
  }

  //
  // Assign a new label to this node.
  //
  void XMLElement::set_label(const std::string& lbl) {
    label_ = lbl;
  }

  //
  // Return the number of attributes held by this element.
  //
  int XMLElement::count_attributes() const {
    return attributes_->size();
  }

  //
  // Return the attribute list held by this element.
  //
  RefPtr<XMLAttributeList> XMLElement::attributes() {
    return attributes_;
  }

  //
  // Assign an attribute list to this element.
  //
  void XMLElement::set_attributes(RefPtr<XMLAttributeList> att) {
    attributes_ = att;
  }

  //
  // Return the number of children held by this element.
  //
  int XMLElement::count_children() const {
    return children_.size();
  }

  //
  // Get the child vertex with the given index.
  //
  RefPtr<XMLElement> XMLElement::get_child(int index) {
    if(index < 0 || index >= int(children_.size()))
      throw MyException("XMLElement::get_child:  Invalid index");
    else
      return children_[index];
  }
  
  //
  // Get the first instance of a child with the given label
  //
  RefPtr<XMLElement> XMLElement::get_child(const std::string& lbl) {
    // Couldn't quite figure out how to easily create a predicate w/o first becoming an
    // expert in the STL, so I programmed a simple search algorithm myself...
    //child_iterator index = std::find_if(children_.begin(), children_.end(), predicate );
    
    bool done = 0;
    int index = 0;
    while(!done && index < int(children_.size())){
      if(strcmp(children_[index]->label().c_str(),lbl.c_str())==0){
        done=1;
      } else {
        index++;
      }
    }
        
    if(!done)
      throw MyException("XMLElement::get_child:  child label not found");
    else
      return children_[index];
  }

  //
  // Add a child vertex to this node.
  //
  void XMLElement::add_child(RefPtr<XMLElement> kid) {
    // first, make sure that this element is not a current child.
    // Note, this checks whether a given pointer value has already been stored,
    // not whether a child with a given label has been stored. Child labels do not
    // have to be unique
    if(std::find(children_.begin(), children_.end(), kid) != children_.end())
      return;

    // make sure the kid is not a NULL pointer.
    if(kid == NULL)
      throw MyException("XMLElement::add_child:  NULL pointer given.");
    // make sure we are not adding this vertex as a kid.
    if(kid == this)
      throw MyException
	("XMLElement::add_child:  Can't add self as child.");

    // and, finally, make sure we don't have a cyclic relationship.
    std::set< RefPtr<XMLElement> > seen;
    recurse(kid, seen);

    // If we get here, it is safe to add the child vertex.
    children_.push_back(kid);
  }

 //
  // Add a non-unique child vertex to this node.
  //
  void XMLElement::add_child_rpt(RefPtr<XMLElement> kid) {
    // The same as add_child, only it allows repeating children (Kh.S. July, 2010)
  

    // make sure the kid is not a NULL pointer.
    if(kid == NULL)
      throw MyException("XMLElement::add_child:  NULL pointer given.");
    // make sure we are not adding this vertex as a kid.
    if(kid == this)
      throw MyException
	("XMLElement::add_child:  Can't add self as child.");

    // and, finally, make sure we don't have a cyclic relationship.
    std::set< RefPtr<XMLElement> > seen;
    recurse(kid, seen);

    // If we get here, it is safe to add the child vertex.
    children_.push_back(kid);
  }


  //
  // Clear all children from this element.
  //
  void XMLElement::clear_children() {
    children_.clear();
  }

  //
  // Return the number of lines of content held by this element.
  //
  int XMLElement::count_content() const {
    return content_.size();
  }

  //
  // Return the content line with the given index.
  //
  const std::string& XMLElement::get_content_line(int index) {
    if(index < 0 || index >= int(content_.size()))
      throw MyException("XMLElement::get_content_line:  Invalid index.");
    else
      return content_[index];
  }
  
  //
  // Add a line of content.
  //
  void XMLElement::add_content_line(const std::string& text) {
    content_.push_back(text);
  }

  //
  // Clear all the content lines.
  //
  void XMLElement::clear_content() {
    content_.clear();
  }

  //
  // Private routine to ensure that we don't have a cyclic relationship.
  //
  void XMLElement::recurse(RefPtr<XMLElement> kid,
			      std::set< RefPtr<XMLElement> > seen) {
    // Only step on each node once.
    if(seen.find(kid) != seen.end())
      return;

    // If this element has been seen before, we have a cyclic graph.
    // The same holds true if this element is the same as the current element.
    if(kid == this)
      throw MyException(std::string("XMLElement::add_child:  Adding ") +
			   "this vertex would yield a cyclic relationship");

    // add this kid to the set
    seen.insert(kid);

    int offspring = kid->count_children();
    for(int i = 0; i < offspring; ++i)
      recurse(kid->get_child(i), seen);
  }
