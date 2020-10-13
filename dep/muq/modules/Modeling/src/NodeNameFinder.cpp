#include "MUQ/Modeling/NodeNameFinder.h"

using namespace muq::Modeling;

NodeNameFinder::NodeNameFinder(std::string const& name, Graph const& graph) : name(name), graph(graph) {}

bool NodeNameFinder::operator()(boost::graph_traits<Graph>::vertex_descriptor vertex) const {
  // check the names
  return name.compare(graph[vertex]->name)==0;
}
