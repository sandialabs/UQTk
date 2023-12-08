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
#include "XMLUtils.h"

//
// Dump out XML tree to a file or output stream
//
void dump_xml_tree(RefPtr<XMLElement> tree, const std::string& indentation,std::ostream& outfile) {
  // Print the 'header' for this node.
  outfile << indentation << "<" << tree->label();

  // Dump all attributes associated with this node.
  RefPtr<XMLAttributeList> att = tree->attributes();
  XMLAttributeList::iterator it, end = att->end();
  for(it = att->begin(); it != end; ++it)
    outfile << " " << it->first << "=\"" << it->second << "\"";

  if(tree->count_children() == 0) {
    // If this node has no children, we're done.
    outfile << " />\n";
  }
  else {
    // ... otherwise, we recursively dump all the child nodes.
    outfile << ">\n";
    int children = tree->count_children();
    std::string new_indent = indentation + "  ";
    for(int kid = 0; kid < children; ++kid)
      dump_xml_tree(tree->get_child(kid), new_indent, outfile);
    // .. and print the closing bracket for this element.
    outfile << indentation << "</" << tree->label() << ">\n";
  }
}

