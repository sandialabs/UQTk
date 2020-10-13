#ifndef NODENAMEFINDER_H_
#define NODENAMEFINDER_H_

#include <boost/graph/adjacency_list.hpp>

#include "MUQ/Modeling/WorkGraphNode.h"
#include "MUQ/Modeling/WorkGraphEdge.h"

namespace muq {
  namespace Modeling {

    /// Define a directed graph type
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, std::shared_ptr<WorkGraphNode>, std::shared_ptr<WorkGraphEdge> > Graph;

    /// A helper struct that determines if a node in the graph has a given name
    struct NodeNameFinder {
    public:

      /**
	 @param[in] name We are looking for nodes with this name
	 @param[in] graph A pointer to the graph that stores the nodes
      */
      NodeNameFinder(std::string const& name, Graph const& graph);

      /// Does a given vertex have the same name as the given name
      /**
	 param[in] vertex The vertex of the graph
	 \return true if the names are the same, false if not
      */
      bool operator()(boost::graph_traits<Graph>::vertex_descriptor vertex) const;

      /// We are looking for vertices with this name
      const std::string name;

      /// This graph stores the vertices
      Graph const& graph;
    };
  } // namespace Modeling
} // namespace muq

#endif
