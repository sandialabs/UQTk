/* =====================================================================================
                     The UQ Toolkit (UQTk) version @UQTKVERSION@
                     Copyright (@UQTKYEAR@) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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

