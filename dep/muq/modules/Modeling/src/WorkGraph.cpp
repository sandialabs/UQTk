#include "MUQ/Modeling/WorkGraph.h"

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ConstantPiece.h"

#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Modeling/LinearAlgebra/IdentityOperator.h"

#include <fstream>
#include <algorithm>

#include <boost/algorithm/string.hpp>

// boost graph library includes
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>

using namespace muq::Modeling;

/// A helper struct that determines if an edge has the same input number
struct SameInputDim {
  /**
     @param[in] dim The input number of the edge
     @param[in] graph A pointer to the graph that stores the edges
   */
  SameInputDim(unsigned int const dim, Graph const& graph) : dim(dim), graph(graph) {}

  /**
     @param[in] edge The edge we want to compate the input number to
     \return true if the input number is the same, false otherwise
  */
  bool operator()(boost::graph_traits<Graph>::edge_descriptor edge) {
    return graph[edge]->inputDim == dim;
  }

  /// The input number of the edge
  int dim;

  /// This graph stores the edges
  Graph const& graph;
};

WorkGraph::WorkGraph() {}

WorkGraph::~WorkGraph() {}

std::shared_ptr<WorkGraph> WorkGraph::Clone() const
{
  std::shared_ptr<WorkGraph> newGraph = std::make_shared<WorkGraph>();
  copy_graph(graph, newGraph->graph);
  return newGraph;
}

void WorkGraph::RemoveNode(std::string const& nodeName)
{
  auto nodeDesc = GetNodeIterator(nodeName);

  clear_vertex(*nodeDesc, graph);  //remove edges to the node
  remove_vertex(*nodeDesc, graph); //then remove the node
}

std::vector<std::pair<int,int>> WorkGraph::GetEdges(std::string const& srcName, std::string const& tgtName)
{
  auto nodeDesc = GetNodeIterator(srcName);

  std::vector<std::pair<int,int>> edges;

  boost::graph_traits<Graph>::out_edge_iterator e, e_end;
  for( boost::tie(e, e_end)=out_edges(*nodeDesc, graph); e!=e_end; ++e ) {
    auto vTarget = target(*e, graph);
    if(tgtName == graph[vTarget]->name)
      edges.push_back( std::make_pair( graph[*e]->outputDim, graph[*e]->inputDim ) );
  }

  return edges;
}


unsigned int WorkGraph::NumNodes() const {

  // return the number of vertices
  return boost::num_vertices(graph);
}

unsigned int WorkGraph::NumEdges() const {
  // return the number of edges
  return boost::num_edges(graph);
}


std::vector<std::string> WorkGraph::GetParents(std::string const& name) const
{
  auto v = GetNodeIterator(name);

  std::vector<std::pair<int, std::string> > temp;

  boost::graph_traits<Graph>::in_edge_iterator e, e_end;
  for( boost::tie(e, e_end)=in_edges(*v, graph); e!=e_end; ++e ) {
    auto vSource = source(*e, graph);
    temp.push_back(std::make_pair(graph[*e]->outputDim, graph[vSource]->name));
  }

  // Sort by index
  std::sort(temp.begin(), temp.end());

  std::vector<std::string> output(temp.size());
  for(int i=0; i<temp.size(); ++i)
    output.at(i) = temp.at(i).second;

  return output;
}

std::string WorkGraph::GetParent(std::string const& name, int inputIndex) const
{
  auto v = GetNodeIterator(name);

  boost::graph_traits<Graph>::in_edge_iterator e, e_end;
  for( boost::tie(e, e_end)=in_edges(*v, graph); e!=e_end; ++e ) {
    if(graph[*e]->inputDim == inputIndex){
        auto vSource = source(*e, graph);
        return graph[vSource]->name;
    }
  }

  return "";
}



std::vector<std::string> WorkGraph::GetChildren(std::string const& name) const
{
  auto v = GetNodeIterator(name);

  std::vector<std::string> output;

  boost::graph_traits<Graph>::out_edge_iterator e, e_end;
  for( boost::tie(e, e_end)=out_edges(*v, graph); e!=e_end; ++e ) {
    auto vSource = target(*e, graph);
    output.push_back(graph[vSource]->name);
  }

  return output;
}


std::vector<std::pair<std::string, int> > WorkGraph::GetInputNames() const
{
  auto inputs = GraphInputs();

  std::vector<std::pair<std::string,int>> names(inputs.size());
  for(int i=0; i<inputs.size(); ++i)
    names.at(i) = std::make_pair(graph[inputs.at(i).first]->name, inputs.at(i).second);

  return names;
}

std::vector<std::pair<std::string, int> > WorkGraph::GetOutputNames() const
{
  auto outputs = GraphOutputs();

  std::vector<std::pair<std::string,int>> names(outputs.size());
  for(int i=0; i<outputs.size(); ++i)
    names.at(i) = std::make_pair(graph[outputs.at(i).first]->name, outputs.at(i).second);

  return names;
}



bool WorkGraph::HasNode(std::string const& name) const {

  // create a node iterator
  boost::graph_traits<Graph>::vertex_iterator iter;

  // check if we have the node
  return HasNode(iter, name);
}

bool WorkGraph::HasNode(boost::graph_traits<Graph>::vertex_iterator& iter, std::string const& name) const {

  // try to find the node with this name
  iter = GetNodeIterator(name);

  // the node exists if the iterator is not the end
  return iter!=vertices(graph).second;
}

void WorkGraph::AddNode(std::shared_ptr<WorkPiece> input, std::string const& name) {

  if(HasNode(name))
    throw std::logic_error("Could not add node \"" + name + "\" to graph.  A node with that name already exists." );

  // add a node to the graph
  auto node = add_vertex(graph);

  graph[node] = std::make_shared<WorkGraphNode>(input, name);
}

void WorkGraph::AddEdge(std::string const& nameFrom, unsigned int const outputDim, std::string const& nameTo, unsigned int const inputDim) {
  // get iterators to the upstream and downstream nodes (make sure they exist)
  auto itFrom = GetNodeIterator(nameFrom);
  if(itFrom==vertices(graph).second)
    throw std::logic_error("Could not add an edge from \"" + nameFrom + "\" to \"" + nameTo + "\" because the source node \"" + nameFrom + "\" does not exist in the graph.");

  auto itTo = GetNodeIterator(nameTo);
  if(itTo==vertices(graph).second)
    throw std::logic_error("Could not add an edge from \"" + nameFrom + "\" to \"" + nameTo + "\" because the target node \"" + nameTo + "\" does not exist in the graph.");

  // the number of inputs and outputs
  const int numOutputs = graph[*itFrom]->piece->numOutputs;
  const int numInputs = graph[*itTo]->piece->numInputs;

  // either we don't know the number of outputs from "nameFrom" or the output dimension is less than the number of outputs
  if( numOutputs>=0 && outputDim>=numOutputs )
    throw std::logic_error("Could not add an edge from output " + std::to_string(outputDim) + "\" of \"" + nameFrom + "\" to input " + std::to_string(inputDim) + " of \"" + nameTo + "\" because node \"" + nameFrom + "\" only has " + std::to_string(numOutputs) + " outputs.");

  // either we don't know the number of inputs to "nameTo" or the input dimension is less than the number of inputs
  if( numInputs>=0 && inputDim>=numInputs )
    throw std::logic_error("Could not add an edge from output " + std::to_string(outputDim) + "\" of \"" + nameFrom + "\" to input " + std::to_string(inputDim) + " of \"" + nameTo + "\" because node \"" + nameTo + "\" only has " + std::to_string(numInputs) + " inputs.");

  // the input/output type
  const std::string inType = graph[*itTo]->piece->InputType(inputDim);
  const std::string outType = graph[*itFrom]->piece->OutputType(outputDim);

  // either we don't know the input and/or output type or they match
  if(inType.compare("")!=0 && // we don't know the input type
     outType.compare("")!=0 && // we don't know the output type
     inType.compare(outType)!=0 ) { // the types must match
    std::cerr << std::endl << "ERROR: Types do not match in 'WorkGraph::AddEdge'.  The input type node '" << nameTo << "' is " << graph[*itTo]->piece->InputType(inputDim) << " but the output type for node '" << nameFrom << "' is " << graph[*itFrom]->piece->OutputType(outputDim) << std::endl << std::endl;
    assert(inType.compare(outType)==0); // the types must match
  }


  // Check to see if the nodes are ModPieces and then check the sizes
  auto modPieceTo = std::dynamic_pointer_cast<ModPiece>(graph[*itTo]->piece);
  auto modPieceFrom = std::dynamic_pointer_cast<ModPiece>(graph[*itFrom]->piece);
  if((modPieceTo) && (modPieceFrom)){
    if(modPieceFrom->outputSizes(outputDim) != modPieceTo->inputSizes(inputDim))
      throw std::logic_error("Could not add an edge from output " + std::to_string(outputDim) + "\" of \"" + nameFrom + "\" to input " + std::to_string(inputDim) + " of \"" + nameTo + "\".  The output of \"" + nameFrom + "\" has size " + std::to_string(modPieceFrom->outputSizes(outputDim)) + " but the input of \"" + nameTo + "\" has size " + std::to_string(modPieceTo->inputSizes(inputDim)) + ".");
  }
  // remove any other edge going into dimension inputDim of the nameTo node
  boost::remove_in_edge_if(*itTo, SameInputDim(inputDim, graph), graph);

  // try to add the new edge, if an edge already exists, notFound will be false and we need to delete the current edge first
  auto temp = boost::add_edge(*itFrom, *itTo, graph);

  // set the edge to have the current dimension
  graph[temp.first] = std::make_shared<WorkGraphEdge>(outputDim, inputDim);
}

boost::graph_traits<Graph>::vertex_iterator WorkGraph::GetNodeIterator(std::string const& name) const {

  // get iterators to the begining and end of the graph
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  boost::tie(v, v_end) = vertices(graph);

  // return the iterator with this name (it is end if that node does not exist)
  auto res = std::find_if(v, v_end, NodeNameFinder(name, graph));

  return res;
}
boost::graph_traits<Graph>::vertex_iterator WorkGraph::GetNodeIterator(std::shared_ptr<WorkPiece> piece) const {

  // get iterators to the begining and end of the graph
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  boost::tie(v, v_end) = vertices(graph);

  // return the iterator with this name (it is end if that node does not exist)
  auto res = std::find_if(v, v_end, [this,piece](boost::graph_traits<Graph>::vertex_descriptor vertex)->bool { return piece==this->graph[vertex]->piece; } );

  return res;

}

std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > WorkGraph::GraphOutputs() const {
  // create an empty vector to hold outputs
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > outputs;

    // loop through the vertices
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  for( std::tie(v, v_end)=vertices(graph); v!=v_end; ++v ) {
    // a vector of the outputs that are set
    std::vector<int> isSet;

    // number of outputs (negative indcates we don't know)
    const int numOutputs = graph[*v]->piece->numOutputs;

    // if possible, reserve memory for the outputs that are set (the size reserved is the total number of outputs, however, if it is negative we don't know how many outputs there are so we reserve whatever the compiler did by default...)
    isSet.reserve(std::max((int)isSet.capacity(), numOutputs));

    // for each vertex, loop over the input nodes and figure out if the outputs are set
    boost::graph_traits<Graph>::out_edge_iterator e, e_end;
    for( tie(e, e_end)=out_edges(*v, graph); e!=e_end; ++e ) {
      // we have an input, so store it in the vector
      isSet.push_back(graph[*e]->outputDim);
    }

    // if an input to this ModPiece is not set, it will be stored as an output to the graph
    for( int i=0; i<numOutputs; ++i ) {
      if( std::find(isSet.begin(), isSet.end(), i)==isSet.end() ) { // if the output is not set ..
	      // ... store this vertex and the output number
	      outputs.push_back(std::make_pair(*v, i));
      }
    }

    // if we don't know the number of outputs
    if( numOutputs<0 ) {
      outputs.push_back(std::make_pair(*v, -1));

      if( isSet.size()>0 ) { // if some outputs have been set ...
	      // the maximum output number that is set
	      const unsigned int maxOut = *std::max_element(isSet.begin(), isSet.end());

	      // loop through all the outputs less than the max input
	      for( unsigned int i=0; i<maxOut; ++i ) {
	        if( std::find(isSet.begin(), isSet.end(), i)==isSet.end() ) { // if an output is not set ...
	          // ... it must be a graph output
	          outputs.push_back(std::make_pair(*v, i));
	        }
	      }
      }
    }
  }

  return outputs;
}

std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > WorkGraph::GraphInputs() const {
  // create an empty vector to hold inputs
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > inputs;

  // loop through the vertices
  boost::graph_traits<Graph>::vertex_iterator v, v_end;
  for( std::tie(v, v_end)=vertices(graph); v!=v_end; ++v ) {
    // a vector of the inputs that are set
    std::vector<int> isSet;

    // number of inputs (negative indcates we don't know)
    const int numInputs = graph[*v]->piece->numInputs;

    // if possible, reserve memory for the inputs that are set (the size reserved is the total number of inputs, however, if it is negative we don't know how many inputs there are so we reserve whatever the compiler did by default...)
    isSet.reserve(std::max((int)isSet.capacity(), numInputs));

    // for each vertex, loop over the input nodes and figure out if the inputs are set
    boost::graph_traits<Graph>::in_edge_iterator e, e_end;
    for( tie(e, e_end)=in_edges(*v, graph); e!=e_end; ++e ) {
      // we have an input, so store it in the vector
      isSet.push_back(graph[*e]->inputDim);
    }

    // if an input to this ModPiece is not set, it will be stored as an input to the graph
    for( int i=0; i<numInputs; ++i ) {
      if( std::find(std::begin(isSet), std::end(isSet), i)==isSet.end() ) { // if the input is not set ..
	      // ... store this vertex and the input number
	      inputs.push_back(std::make_pair(*v, i));
      }
    }

    // if we don't know the number of inputs
    if( numInputs<0 ) {
      inputs.push_back(std::make_pair(*v, -1));

      if( isSet.size()>0 ) { // if some inputs have been set ...
	      // the maximum input number that is set
	      const unsigned int maxIn = *std::max_element(isSet.begin(), isSet.end());

	      // loop through all the inputs less than the max input
	      for( unsigned int i=0; i<maxIn; ++i ) {
	        if( std::find(isSet.begin(), isSet.end(), i)==isSet.end() ) { // if an input is not set ...
	          // ... it must be a graph input
	          inputs.push_back(std::make_pair(*v, i));
	        }
	      }
      }
    }
  }

  return inputs;
}

bool WorkGraph::HasEdge(boost::graph_traits<Graph>::vertex_descriptor const& vOut, boost::graph_traits<Graph>::vertex_descriptor const& vIn, int const inputDim) const {
  // iteraters through the output edges of the input node
  boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
  boost::tie(ei, ei_end) = boost::out_edges(vOut, graph);

  for( ; ei!=ei_end ; ++ei ) {  // loop through the output edges
    if( target(*ei, graph)==vIn && graph[*ei]->inputDim==inputDim ) { // if the target is the input node and the input dimension is the same
      return true;
    }
  }

  return false;
}

void WorkGraph::RecursiveCut(const boost::graph_traits<Graph>::vertex_descriptor& vOld, const boost::graph_traits<Graph>::vertex_descriptor& vNew, std::shared_ptr<WorkGraph>& newGraph) const {
  // a map from the source ID to a pair: <source vertex, vector of edges from that vertex to this one>
  std::map<unsigned int, std::pair<boost::graph_traits<Graph>::vertex_descriptor, std::vector<boost::graph_traits<Graph>::in_edge_iterator> > > sources;

  // the input edges into vOld
  boost::graph_traits<Graph>::in_edge_iterator e, e_end;

  for( tie(e, e_end)=in_edges(vOld, graph); e!=e_end; ++e ) { // loop through the input edges
    auto v = source(*e, graph);
    const unsigned int id = graph[v]->piece->ID();

    auto it = sources.find(id);
    if( it==sources.end() ){
      sources[id] = std::pair<boost::graph_traits<Graph>::vertex_descriptor, std::vector<boost::graph_traits<Graph>::in_edge_iterator> >(v, std::vector<boost::graph_traits<Graph>::in_edge_iterator>(1, e));
    } else {
      sources[id].second.push_back(e);
    }
  }

  // loop through the sources
  for( auto it : sources ) {
    // the upstream node iterator
    auto v = it.second.first;

    if( Constant(v) && std::dynamic_pointer_cast<ConstantPiece>(graph[v]->piece)==nullptr ) { // if this node is constant but is not already a muq::Modeling::ConstantPiece
      // get the output values for this node
      const std::vector<boost::any>& outputs = GetConstantOutputs(v);

      // create a ConstantPiece node for this input
      auto nextV = boost::add_vertex(newGraph->graph);
      newGraph->graph[nextV] = std::make_shared<WorkGraphNode>(std::make_shared<ConstantPiece>(outputs), graph[v]->name+"_fixed");

      // loop through the edges from the source to this node
      for( auto e : it.second.second ) {
	if( !newGraph->HasEdge(nextV, vNew, graph[*e]->inputDim) ) { // if edge does not exist ...
	  // ... add the edge from this node to the existing node
	  auto nextE = boost::add_edge(nextV, vNew, newGraph->graph);
	  newGraph->graph[nextE.first] = std::make_shared<WorkGraphEdge>(graph[*e]->outputDim, graph[*e]->inputDim);
	}
      }

      // move to the next source
      continue;
    }

    // loop through the edges from the (nonconstant) source to this node
    for( auto e : it.second.second ) {
      boost::graph_traits<Graph>::vertex_descriptor nextV;
      boost::graph_traits<Graph>::vertex_iterator ind;

      if( newGraph->HasNode(ind, graph[v]->name) ) { // if the node already exists in the new graph ...
	      // ... the next node is that node
	      nextV = *ind;
      } else { // if not ...
	      // ... copy the source node over to the new graph and make an edge
	      nextV = boost::add_vertex(newGraph->graph);
	      newGraph->graph[nextV] = graph[v];
      }

      if( !newGraph->HasEdge(nextV, vNew, graph[*e]->inputDim ) ) { // if the edge does not already exist ...
	      // add the edge from this node to the existing node,
	      auto nextE = boost::add_edge(nextV, vNew, newGraph->graph);
	      newGraph->graph[nextE.first] = std::make_shared<WorkGraphEdge>(graph[*e]->outputDim, graph[*e]->inputDim);
      }

      // recurse down a step further
      RecursiveCut(v, nextV, newGraph);
    }
  }
}

std::shared_ptr<WorkGraph> WorkGraph::DependentCut(std::string const& nameOut) const {
  // create a new graph
  auto newGraph = std::make_shared<WorkGraph>();

  // the an iterator to the output node
  auto oldV = GetNodeIterator(nameOut);

  // if the desired node is constant
  if( Constant(nameOut) ) {
    // get the output values for this node
    const std::vector<boost::any>& outputs = GetConstantOutputs(nameOut);

    // create a ConstantPiece node for this input
    auto nextV = boost::add_vertex(newGraph->graph);
    newGraph->graph[nextV] = std::make_shared<WorkGraphNode>(std::make_shared<ConstantPiece>(outputs), graph[*oldV]->name+"_fixed");

    // return a graph with only one (constant node)
    return newGraph;
  }

  // add a new node to the output graph
  auto newV = boost::add_vertex(newGraph->graph);

  // copy the output node
  newGraph->graph[newV] = graph[*oldV];

  // recurse through graph
  RecursiveCut(*oldV, newV, newGraph);

  return newGraph;
}

std::shared_ptr<WorkGraphPiece> WorkGraph::CreateWorkPiece(std::string const& node) const {

      // make sure we have the node
  assert(HasNode(node));

  // trime the extraneous branches from the graph
  auto newGraph = DependentCut(node);
  assert(newGraph);

  // get the inputs to the cut graph
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > inputs = newGraph->GraphInputs();

  // the name of each input
  std::vector<std::string> inputNames;
  inputNames.reserve(inputs.size()); // reserve the size

  // loop through the input nodes and create each input name
  for( auto it=inputs.begin(); it!=inputs.end(); ++it ) {
    // make sure the input number is known
    if( newGraph->graph[it->first]->piece->numInputs<0 ) {
        std::cerr << std::endl << "ERROR: Cannot create WorkGraphPiece if one of the nodes has an unknown number of inputs.  Node \"" << newGraph->graph[it->first]->name<< "\" does not specify the number of inputs. " << std::endl << std::endl;

      assert(newGraph->graph[it->first]->piece->numInputs>=0);
    }

    std::stringstream temp;
    temp << newGraph->graph[it->first]->name << "_";
    temp << it->second;

    inputNames.push_back(temp.str());
  }

  // the constant pieces that will hold the inputs
  std::vector<std::shared_ptr<ConstantPiece> > constantPieces(inputs.size());

  // the input types
  std::map<unsigned int, std::string> inTypes;

  assert(inputNames.size()==inputs.size());
  assert(constantPieces.size()==inputs.size());
  for( unsigned int i=0; i<inputs.size(); ++i ) { // loop over each input
    const std::string inType = newGraph->graph[inputs.at(i).first]->piece->InputType(inputs.at(i).second, false);
    if( inType.compare("")!=0 ) {
      inTypes[i] = inType;
    }

    // create a constant WorkPiece to hold the input (it is empty for now) and add it to the new graph
    constantPieces[i] = std::make_shared<ConstantPiece>();
    newGraph->AddNode(constantPieces[i], inputNames[i]);
    newGraph->AddEdge(inputNames[i], 0, newGraph->graph[inputs[i].first]->name, inputs[i].second);
  }

  // Look for the original node name
  auto outNode = newGraph->GetNodeIterator(node);

  // If we didn't find the original node, look for the fixed one
  if(outNode == vertices(newGraph->graph).second){
      std::string node_fixed = node + "_fixed";
      outNode = newGraph->GetNodeIterator(node_fixed);
      assert(outNode != vertices(newGraph->graph).second);
  }

  //return std::make_shared<WorkGraphPiece>(newGraph->graph, constantPieces, inputName, inTypes, newGraph->graph[*outNode]->piece);
  return std::make_shared<WorkGraphPiece>(newGraph, constantPieces, inputNames, inTypes, newGraph->graph[*outNode]->piece);

}


std::shared_ptr<ModGraphPiece> WorkGraph::CreateModPiece(std::string const& node, std::vector<std::string> const& inNames) const {

  // make sure we have the node
  assert(HasNode(node));

  // trime the extraneous branches from the graph
  auto newGraph = DependentCut(node);
  assert(newGraph);

  // get the inputs to the cut graph
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > inputs;
  if( inNames.size()>0 ) {
    for( unsigned int i=0; i<inNames.size(); ++i ) {
      boost::graph_traits<Graph>::vertex_iterator it;
      const bool has = newGraph->HasNode(it, inNames[i]);
      assert(has);

      // get the number of inputs
      const int numInputs = newGraph->graph[*it]->piece->numInputs;
      std::vector<int> isSet;
      isSet.reserve(numInputs);
      // for each vertex, loop over the input nodes and figure out if the inputs are set
      boost::graph_traits<Graph>::in_edge_iterator e, e_end;
      for( tie(e, e_end)=in_edges(*it, newGraph->graph); e!=e_end; ++e ) {
        // we have an input, so store it in the vector
        isSet.push_back(graph[*e]->inputDim);
      }

      // for each vertex, loop over the input nodes and figure out if the inputs are set
      for( int j=0; j<numInputs; ++j ) {
        if( std::find(std::begin(isSet), std::end(isSet), j)==isSet.end() ) { // if the input is not set ..
          inputs.push_back(std::pair<boost::graph_traits<Graph>::vertex_descriptor, int>(*it, j));
        }
      }
    }
  } else {
    inputs = newGraph->GraphInputs();
  }

  // the name of each input
  std::vector<std::string> inputNames;
  inputNames.reserve(inputs.size()); // reserve the size

  // loop through the input nodes and create each input name
  for( auto it=inputs.begin(); it!=inputs.end(); ++it ) {
    // make sure the input number is known
    if( newGraph->graph[it->first]->piece->numInputs<0 ) {
        std::cerr << std::endl << "ERROR: Cannot create WorkGraphPiece if one of the nodes has an unknown number of inputs.  Node \"" << newGraph->graph[it->first]->name<< "\" does not specify the number of inputs. " << std::endl << std::endl;

      assert(newGraph->graph[it->first]->piece->numInputs>=0);
    }

    std::stringstream temp;
    temp << newGraph->graph[it->first]->name << "_";
    temp << it->second;

    inputNames.push_back(temp.str());
  }

  // the constant pieces that will hold the inputs
  std::vector<std::shared_ptr<ConstantVector> > constantPieces(inputs.size());

  assert(inputNames.size()==inputs.size());
  assert(constantPieces.size()==inputs.size());

  for( unsigned int i=0; i<inputs.size(); ++i ) { // loop over each input
    // create a constant WorkPiece to hold the input (it is empty for now) and add it to the new graph
    auto modIn = std::dynamic_pointer_cast<ModPiece>(newGraph->graph[inputs.at(i).first]->piece);
    assert(modIn);

    constantPieces.at(i) = std::make_shared<ConstantVector>(Eigen::VectorXd::Zero(modIn->inputSizes(inputs.at(i).second)));
    newGraph->AddNode(constantPieces.at(i), inputNames.at(i));
    newGraph->AddEdge(inputNames.at(i), 0,
                      newGraph->graph[inputs.at(i).first]->name, inputs.at(i).second);
  }

  // Look for the original node name
  auto outNode = newGraph->GetNodeIterator(node);

  // If we didn't find the original node, look for the fixed one
  if(outNode == vertices(newGraph->graph).second){
      std::string node_fixed = node + "_fixed";
      outNode = newGraph->GetNodeIterator(node_fixed);
      assert(outNode != vertices(newGraph->graph).second);
  }

  //return std::make_shared<WorkGraphPiece>(newGraph->graph, constantPieces, inputName, inTypes, newGraph->graph[*outNode]->piece);
  return std::make_shared<ModGraphPiece>(newGraph, constantPieces, inputNames, std::dynamic_pointer_cast<ModPiece>(newGraph->graph[*outNode]->piece));
}

/** Print the nodes and edges of this graph to std::cout. */
void WorkGraph::Print(std::ostream& fout) const
{
  // first, list the vertices
  fout << "\nNodes:\n";
  boost::graph_traits<Graph>::vertex_iterator v, v_end;

  for (std::tie(v, v_end) = vertices(graph); v != v_end; v++) {
    fout << "\t" << graph[*v]->name << std::endl;
  }

  // now, list the edges
  fout << "Edges:\n";
  boost::graph_traits<Graph>::edge_iterator e, e_end;
  for (std::tie(e, e_end) = edges(graph); e != e_end; e++) {
      fout << "\t" <<
      graph[source(*e,graph)]->name << "[" << graph[*e]->outputDim << "] -> " <<
      graph[target(*e,graph)]->name << "[" << graph[*e]->inputDim << "]\n";
  }
  fout << "\n";
}



std::shared_ptr<WorkPiece> WorkGraph::GetPiece(std::string const& name)
{
  return graph[*GetNodeIterator(name)]->piece;
}
std::shared_ptr<WorkPiece> WorkGraph::GetPiece(boost::graph_traits<Graph>::vertex_descriptor it)
{
  return graph[it]->piece;
}

struct TrueOp {
  template<typename T>
  bool operator()(const T& in)
  {
    return true;
  }
};

std::string WorkGraph::GetName(std::shared_ptr<WorkPiece> piece) const
{
  return graph[*GetNodeIterator(piece)]->name;
}

void WorkGraph::BindNode(std::string const& nodeName, std::vector<boost::any> const& x)
{
  // find the node
  auto nodeDesc = GetNodeIterator(nodeName);

  // next, delete all incoming edges
  boost::remove_in_edge_if(*nodeDesc, TrueOp(), graph);

  // try to cast the node as a ModPiece
  auto mod = std::dynamic_pointer_cast<ModPiece>((graph)[*nodeDesc]->piece);
  if( mod ) {
    // replace the ModPiece ptr
    std::vector<Eigen::VectorXd> vec(x.size());
    for( unsigned int i=0; i<x.size(); ++i ) {
      vec[i] = boost::any_cast<Eigen::VectorXd const>(x[i]);
    }
    (graph)[*nodeDesc]->piece = std::make_shared<ConstantVector>(vec);
  } else {
    // replace the WorkPiece ptr
    (graph)[*nodeDesc]->piece = std::make_shared<ConstantPiece>(x);
  }
}

void WorkGraph::BindEdge(std::string const& nodeName,
                         unsigned int       inputDim,
                         boost::any  const& x)
{
  auto nodeDesc = GetNodeIterator(nodeName);

  // iterate through the input edges to determine if this edge already exists, remove it if it does
  boost::graph_traits<Graph>::in_edge_iterator e, e_end;

  for (std::tie(e, e_end) = boost::in_edges(*nodeDesc, graph); e != e_end; e++) {
    if ((graph)[*e]->inputDim == inputDim) {
      boost::remove_edge(*e, graph);
      break;
    }
  }

  std::string newName = nodeName + "_FixedInput" + std::to_string(inputDim);
  auto newPiece = std::make_shared<ConstantPiece>(std::vector<boost::any>(1,x));
  AddNode(newPiece, newName);
  AddEdge(newName, 0, nodeName, inputDim);
}

class MyVertexWriter {
public:

  MyVertexWriter(Graph const& graph) : graph(graph) {}

  void operator()(std::ostream& out, const boost::graph_traits<Graph>::vertex_descriptor& v) const {
    int status;

    // the WorkPiece associated with this node
    auto workPtr = graph[v]->piece;

    // the name of this work piece
    const std::string nodeName = workPtr->Name();

    // style for the node visualization
    const std::string style = "colorscheme=pastel16,color=2, style=filled";

    // label the node
    out << "[label=\"" << graph[v]->name << " : " << nodeName << "\", " << style << "]";
  }

private:
  // the graph we are visualizing
  Graph const& graph;
};

class MyEdgeWriter {
public:

  MyEdgeWriter(Graph const& graph) : graph(graph) {}

  void operator()(std::ostream& out, const boost::graph_traits<Graph>::edge_descriptor& e) const {
    const unsigned int inputDim = graph[e]->inputDim;

    const unsigned int outputDim = graph[e]->outputDim;

    // first, write the name as the label
    out << "[label=\" [out, in]: [" << outputDim << ", " << inputDim << "]\"]";
  }

private:
  // the graph we are visualizing
  Graph const& graph;
};

class MyGraphWriter {
public:

  MyGraphWriter(Graph const& graph) : graph(graph) {}

  void operator()(std::ostream& out) const {
    out << "splines = true;" << std::endl;
  }

private:
  Graph const& graph;
};

void WorkGraph::Visualize(std::string const& filename) const {
  // split the graph (name and extension)
  std::vector<std::string> strs;
  boost::split(strs, filename, boost::is_any_of("."));

  // is the extension something we expect?
  const bool knownExtension = (strs.end()-1)->compare("png") || (strs.end()-1)->compare("jpg") || (strs.end()-1)->compare("tif") || (strs.end()-1)->compare("eps") || (strs.end()-1)->compare("pdf") || (strs.end()-1)->compare("svg");

  // open a file stream
  std::ofstream fout;

  const std::string tempname = *strs.begin() + "_temp.dot";
  if( knownExtension ) { // if we know the file extension ...
    // ... open a temporary file
    fout.open(tempname);
  } else { // otherwise ...
    // ... just open the file as the user asked
    fout.open(filename.c_str());
  }

  typedef std::map<boost::graph_traits<Graph>::vertex_descriptor, size_t> IndexMap;
  IndexMap mapIndex;
  boost::associative_property_map<IndexMap> propmapIndex(mapIndex);

  int vertexNum = 0;
  BGL_FORALL_VERTICES(v, graph, Graph) {
    put(propmapIndex, v, vertexNum++);
  }

  boost::write_graphviz(fout, graph, MyVertexWriter(graph), MyEdgeWriter(graph), MyGraphWriter(graph), propmapIndex);

  fout.seekp(-2, std::ios_base::cur); //back up so the rest is inside the brackets

  // loop over all the inputs and draw them
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > graphInputs = GraphInputs();

  int in = 0;
  for( auto aPair : graphInputs ) { // loop through the graph inputs
    vertexNum++;
    if( aPair.second<0 ) { // if we do not know the input number
      fout << vertexNum << "[label=\"Unfixed input\", shape=invhouse,colorscheme=pastel13,color=1, style=filled];" << std::endl;
      fout << vertexNum << "->" << propmapIndex[aPair.first] << std::endl;
    } else { // if we know the input number
      fout << vertexNum << "[label=\"Input #" << in << "\", shape=invhouse,colorscheme=pastel13,color=1, style=filled];" << std::endl;
      fout << vertexNum << "->" << propmapIndex[aPair.first] << "[label=\" in: " << aPair.second << "\"];" << std::endl;
      ++in;
    }
  }

  // loop over all the outputs and draw them
  std::vector<std::pair<boost::graph_traits<Graph>::vertex_descriptor, int> > graphOutputs = GraphOutputs();

  int out = 0;
  for( auto aPair : graphOutputs ) { // loop through the graph outputs
    vertexNum++;
    if( aPair.second<0 ) { // if we do not know the output number
      fout << vertexNum << "[label=\"Unfixed output\", shape=box,colorscheme=pastel16,color=1, style=filled];" << std::endl;
      fout << propmapIndex[aPair.first]  << "->" << vertexNum << std::endl;
    } else { // if we know the input number
      fout << vertexNum << "[label=\"Output #" << out << "\", shape=box,colorscheme=pastel16,color=1, style=filled];" << std::endl;
      fout << propmapIndex[aPair.first]  << "->" << vertexNum << "[label=\" out: " << aPair.second << "\"];" << std::endl;
      ++out;
    }
  }

  fout << "}" << std::endl;

  // close the file stream
  fout.close();

  if( knownExtension ) {
    // move the *.dot file into the *.[whatever extension file]
    std::system(("dot -T" + *(strs.end() - 1) + " " + tempname + " -o " + filename).c_str());

    // remove the temporary *.dot file
    std::system(("rm "+tempname).c_str());
  }
}

std::vector<boost::any> const& WorkGraph::GetConstantOutputs(std::string const& node) const {
  // make sure the node indeed cosntant
  assert(Constant(node));

  return GetConstantOutputs(*GetNodeIterator(node));
}

std::vector<boost::any>& WorkGraph::GetConstantOutputs(boost::graph_traits<Graph>::vertex_descriptor const& node) const {
  // make sure the node indeed cosntant
  assert(Constant(node));

  // the WorkPiece associated with this node
  auto work = graph[node]->piece;

  // make sure we know the number of inputs
  assert(work->numInputs>=0);

  // if the node has no inputs
  if( work->numInputs==0 ) {
    work->Evaluate();
    return work->outputs;
  }

  // get iterators to the input edges
  boost::graph_traits<Graph>::in_edge_iterator e, e_end;

  // a map from the WorkPiece ID number to a pair: <upstream vertex descriptor, vector of pairs: <ouput supply, input supplied> >
  std::map<unsigned int, std::pair<boost::graph_traits<Graph>::vertex_descriptor, std::vector<std::pair<unsigned int, unsigned int> > > > inMap;

  // loop through the input edges
  for( tie(e, e_end)=in_edges(node, graph); e!=e_end; ++e ) {
    // get the WorkPiece id number, the output that it supplies, and the input that receives it
    const unsigned int id = graph[source(*e, graph)]->piece->ID();
    const unsigned int inNum = graph[*e]->inputDim;
    const unsigned int outNum = graph[*e]->outputDim;

    // add to the list of inputs supplyed by the WorkPiece with this id
    auto it = inMap.find(id);
    if( it==inMap.end() ) {
      inMap[id] = std::pair<boost::graph_traits<Graph>::vertex_descriptor, std::vector<std::pair<unsigned int, unsigned int> > >(boost::source(*e, graph), std::vector<std::pair<unsigned int, unsigned int> >(1, std::pair<unsigned int, unsigned int>(outNum, inNum)));
    } else {
      inMap[id].second.push_back(std::pair<unsigned int, unsigned int>(outNum, inNum));
    }
  }

  // the inputs to this WorkPiece
  boost::any empty(nullptr);
  ref_vector<boost::any> ins(work->numInputs, std::cref(empty));

  // loop through the edges again, now we know which outputs supply which inputs
  for( auto it : inMap ) {
    // get the outputs of the upstream node
    std::vector<boost::any>& upstreamOutputs = GetConstantOutputs(it.second.first);

    // loop through the inputs supplied by the outputs of this upstream node
    for( auto out_in : it.second.second ) {
      // populate the inputs to this WorkPiece
      ins.at(out_in.second) = upstreamOutputs.at(out_in.first);
    }
  }

  // evaluate this node
  work->Evaluate(ins);
  return work->outputs;
}

bool WorkGraph::Constant(std::string const& node) const {
  return Constant(*GetNodeIterator(node));
}

bool WorkGraph::Constant(boost::graph_traits<Graph>::vertex_descriptor const& node) const {
  // the WorkPiece associated with this node
  auto work = graph[node]->piece;

  // if the node has no inputs, it is constant
  if( work->numInputs==0 ) { return true; }

  // if not all of its inputs are set, it is not constant (this also covers the case that we don't know the number of inputs)
  if( in_degree(node, graph)!=work->numInputs ) { return false; }

  // all of the inputs are given by upstream nodes, it is constant if the upstream nodes are constant
  bool constant = true;

  // get iterators to the input edges
  boost::graph_traits<Graph>::in_edge_iterator e, e_end;

  // loop through the input edges
  for( tie(e, e_end)=in_edges(node, graph); e!=e_end; ++e ) {
    // if an upstream node is not constant, then this node is not constant either
    if( !Constant(boost::source(*e, graph)) ) { return false; }
  }

  // all upstream nodes are constant, so this node is alos constant
  return true;
}
