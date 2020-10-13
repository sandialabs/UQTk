#include "MUQ/Modeling/ModGraphPiece.h"

#include "MUQ/Modeling/WorkGraph.h"
#include "MUQ/Utilities/AnyHelpers.h"

#include <boost/range/adaptor/reversed.hpp>
#include <boost/graph/topological_sort.hpp>

#include <Eigen/Core>

using namespace muq::Modeling;
using namespace muq::Utilities;



ModGraphPiece::ModGraphPiece(std::shared_ptr<WorkGraph>                           graph,
                             std::vector<std::shared_ptr<ConstantVector> > const& constantPiecesIn,
                             std::vector<std::string>                      const& inputNames,
                             std::shared_ptr<ModPiece>                            outputPiece) : ModPiece(ConstructInputSizes(constantPiecesIn),
                                                                                                          outputPiece->outputSizes),
                                                                                                 wgraph(graph),
                                                                                                 outputID(outputPiece->ID()),
                                                                                                 constantPieces(constantPiecesIn) {

  // build the run order
  assert(graph);
  boost::topological_sort(wgraph->graph, std::front_inserter(runOrder));

  // each input only needs to loop over its downstream nodes when computing derivatives
  adjointRunOrders.resize(inputSizes.size());
  filtered_graphs.resize(inputSizes.size());

  // compute a run order for each of the inputs so we only have to loop over their downstream nodes
  assert(numInputs==inputNames.size());
  for( unsigned int i=0; i<numInputs; ++i ) { // loop through the inputs

    // get iterators to the begining and end of the graph
    boost::graph_traits<Graph>::vertex_iterator v, v_end;
    boost::tie(v, v_end) = vertices(wgraph->graph);

    // determine the downstream nodes of this input
    DependentPredicate nFilt(*std::find_if(v, v_end, NodeNameFinder(inputNames[i], wgraph->graph)), wgraph->graph);
    DependentEdgePredicate eFilt(nFilt, wgraph->graph);

    // filter the graph, we only care about downstream nodes of this input
    filtered_graphs[i] = std::make_shared<boost::filtered_graph<Graph, DependentEdgePredicate, DependentPredicate> >(wgraph->graph, eFilt, nFilt);

    // build specialized run order for each input dimension
    boost::topological_sort(*filtered_graphs[i], std::back_inserter(adjointRunOrders[i]));
  }
}

Eigen::VectorXi ModGraphPiece::ConstructInputSizes(std::vector<std::shared_ptr<ConstantVector> > const& constantPiecesIn)
{
  Eigen::VectorXi sizes(constantPiecesIn.size());
  for(int i=0; i<constantPiecesIn.size(); ++i){
    sizes(i)  = constantPiecesIn.at(i)->outputSizes(0);
    assert(constantPiecesIn.at(i)->outputSizes.size()==1);
  }
  return sizes;
}

void ModGraphPiece::EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  // set the inputs
  SetInputs(inputs);

  // fill the map from the WorkPiece ID to its outputs
  FillOutputMap();

  // store the result in the output vector
  outputs.resize(valMap[outputID].size());
  for(int i=0; i<outputs.size(); ++i) {
    outputs.at(i) = valMap[outputID].at(i).get();
  }
}

void ModGraphPiece::JacobianImpl(unsigned int                const  wrtOut,
                                 unsigned int                const  wrtIn,
                                 ref_vector<Eigen::VectorXd> const& inputs) {

  // set the inputs
  SetInputs(inputs);

  // fill the map from the WorkPiece ID to its outputs
  FillOutputMap();

  // a map from the WorkPiece ID to a vector holding the cumulative jacobians of that output wrt the specified input
  std::map<unsigned int, std::vector<Eigen::MatrixXd> > jacMap;

  // loop through each downstream node
  for( auto node : boost::adaptors::reverse(adjointRunOrders[wrtIn]) ) {

    std::shared_ptr<ModPiece> nodePiece = std::dynamic_pointer_cast<ModPiece>(filtered_graphs[wrtIn]->operator[](node)->piece);
    assert(nodePiece);

    // the ID of the current node
    const unsigned int nodeID = nodePiece->ID();

    // Initialize the jacobian map for this node
    jacMap[nodeID] = std::vector<Eigen::MatrixXd>(nodePiece->outputSizes.size());

    // get the outputs of this node that impact the specified output node
    const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > & requiredOutNodes = RequiredOutputs(node, wrtIn, wrtOut);

    // remove duplicates
    std::vector<unsigned int> requiredOuts;
    requiredOuts.reserve(requiredOutNodes.size());
    for( auto out : requiredOutNodes ) {
      auto it = std::find(requiredOuts.begin(), requiredOuts.end(), std::get<1>(out));
      if( it==requiredOuts.end() ) {
	       requiredOuts.push_back(std::get<1>(out));
      }
    }

    // Initialize the jacobians that will be stored
    for(int i=0; i<requiredOuts.size(); ++i)
      jacMap[nodeID].at(requiredOuts.at(i)) = Eigen::MatrixXd::Zero(nodePiece->outputSizes(requiredOuts.at(i)), inputSizes(wrtIn));

    // get the inputs for this node --- the input WorkPiece ID, the output number, and the input number
    const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredIns = RequiredInputs(node, wrtIn);

    // if there are no inputs, it is an input node, so set the Jacobian to the identity
    if(requiredIns.size()==0){
      jacMap[nodeID][0] = Eigen::MatrixXd::Identity(nodePiece->outputSizes(0), nodePiece->outputSizes(0));
      assert(jacMap[nodeID].size()==1);
    }else{

      // the inputs to this WorkPiece
      const ref_vector<Eigen::VectorXd>& ins = GetNodeInputs(node);

      // compute the jacobian of each required output wrt each input
      for( auto out : requiredOuts ) {
        // To compute the Jacobian of out, we need to add the add the combination from each input
  	    for( auto in : requiredIns )
  	      jacMap[nodeID][out] += nodePiece->Jacobian(out, std::get<2>(in), ins) * jacMap[std::get<0>(in)][std::get<1>(in)];
      }
    }

  } // loop over run order

  // set the Jacobian for this WorkPiece
  jacobian = jacMap[outputID][wrtOut];
}
//
// void WorkGraphPiece::JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) {
//   // set the inputs
//   SetInputs(inputs);
//
//   // fill the map from the WorkPiece ID to its outputs
//   FillOutputMap();
//
//   // a map from the WorkPiece ID a map from the output number to the action of the jacobian of that output wrt the specified input
//   std::map<unsigned int, std::map<unsigned int, boost::any> > jacActionMap;
//
//   // loop through each downstream node
//   for( auto node : boost::adaptors::reverse(derivRunOrders[wrtIn]) ) {
//     // the ID of the current node
//     const unsigned int nodeID = filtered_graphs[wrtIn]->operator[](node)->piece->ID();
//
//     // get the outputs of this node that depend on the specified input
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredOutNodes = RequiredOutputs(node, wrtIn, wrtOut);
//     // remove duplicates
//     std::vector<unsigned int> requiredOuts;
//     requiredOuts.reserve(requiredOutNodes.size());
//     for( auto out : requiredOutNodes ) {
//       auto it = std::find(requiredOuts.begin(), requiredOuts.end(), std::get<1>(out));
//       if( it==requiredOuts.end() ) {
// 	requiredOuts.push_back(std::get<1>(out));
//       }
//     }
//
//     // get the inputs for this node --- the input WorkPiece ID, the output number, and the input number
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredIns = RequiredInputs(node, wrtIn);
//
//     // the inputs to this WorkPiece
//     const ref_vector<boost::any>& ins = Inputs(node);
//
//     // compute the jacobian of each required output wrt each input
//     for( auto out : requiredOuts ) {
//       if( requiredIns.size()==0 ) {
// 	// if there are no inputs, it is the input so set the Jacobian to the identity
// 	jacActionMap[nodeID][out] = vec;
//       } else {
// 	// initize the jacobian to nothing
// 	jacActionMap[nodeID][out] = boost::none;
//
// 	for( auto in : requiredIns ) {
// 	  // compute the Jacobian with respect to each required input
// 	  graph->operator[](node)->piece->JacobianAction(std::get<2>(in), out, jacActionMap[std::get<0>(in)][std::get<1>(in)], ins);
//
// 	  // use chain rule to get the jacobian wrt to the required input
// 	  jacActionMap[nodeID][out] = algebra->Add(jacActionMap[nodeID][out], *(graph->operator[](node)->piece->jacobianAction));
// 	}
//       }
//     }
//   }
//
//   // set the action of the Jacobian for this WorkPiece
//   jacobianAction = jacActionMap[outputID][wrtOut];
// }
//
// void WorkGraphPiece::JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) {
//     // set the inputs
//   SetInputs(inputs);
//
//   // fill the map from the WorkPiece ID to its outputs
//   FillOutputMap();
//
//   // a map from the WorkPiece ID a map from the output number to the action of the jacobian of that output wrt the specified input
//   std::map<unsigned int, std::map<unsigned int, boost::any> > jacTransActionMap;
//
//   // loop through each downstream node
//   for( auto node : derivRunOrders[wrtIn] ) {
//     // the ID of the current node
//     const unsigned int nodeID = filtered_graphs[wrtIn]->operator[](node)->piece->ID();
//
//     // get the outputs of this node that depend on the specified input
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredOuts = RequiredOutputs(node, wrtIn, wrtOut);
//
//     // get the inputs for this node --- the input WorkPiece ID, the output number, and the input number
//     const std::vector<std::tuple<unsigned int, unsigned int, unsigned int> >& requiredIns = RequiredInputs(node, wrtIn);
//
//     // the inputs to this WorkPiece
//     const ref_vector<boost::any>& ins = Inputs(node);
//
//     for( auto in : requiredIns ) {
//       if( nodeID==outputID ) {
// 	assert(requiredOuts.size()==1);
// 	assert(std::get<1>(requiredOuts[0])==wrtOut);
//
// 	// compute the Jacobian transpose action of the output node
// 	graph->operator[](node)->piece->JacobianTransposeAction(std::get<2>(in), wrtOut, vec, ins);
// 	jacTransActionMap[nodeID][std::get<2>(in)] = *(graph->operator[](node)->piece->jacobianTransposeAction);
//       } else {
// 	// initialize the jacobian transpose action to nothing
// 	jacTransActionMap[nodeID][std::get<2>(in)] = boost::none;
//
// 	// loop through the outputs
// 	for( auto out : requiredOuts ) {
// 	  // compute the jacobian transpose action for this output
// 	  graph->operator[](node)->piece->JacobianTransposeAction(std::get<2>(in), std::get<1>(out), jacTransActionMap[std::get<0>(out)][std::get<2>(out)], ins);
// 	  // add it (chain rule)
// 	  jacTransActionMap[nodeID][std::get<2>(in)] = algebra->Add(jacTransActionMap[nodeID][std::get<2>(in)], *(graph->operator[](node)->piece->jacobianTransposeAction));
// 	}
//       }
//     }
//
//     // if this is the input node
//     if( requiredIns.size()==0 ) {
//       // loop though the outputs
//       for( auto out : requiredOuts ) {
// 	if( jacobianTransposeAction ) { // if the jacobian transpose action has not be initilized ...
// 	  // it is equal to the action of the output
// 	  *jacobianTransposeAction = algebra->Add(*jacobianTransposeAction, jacTransActionMap[std::get<0>(out)][std::get<1>(out)]);
// 	} else {
// 	  // add it to the existing jacobian transpose action (chain rule)
// 	  jacobianTransposeAction = jacTransActionMap[std::get<0>(out)][std::get<1>(out)];
// 	}
//       }
//     }
//   }
// }

void ModGraphPiece::SetInputs(ref_vector<Eigen::VectorXd> const& inputs) {
  // get the inputs and set them to the ConstantPiece nodes
  assert(inputs.size()==constantPieces.size());
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    constantPieces[i]->SetValue(inputs.at(i));
  }
}

std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > ModGraphPiece::InputNodes(boost::graph_traits<Graph>::vertex_descriptor const& node) const {
  // the map of input nodes
  std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > inMap;

  // loop though the input nodes
  boost::graph_traits<Graph>::in_edge_iterator e, e_end;
  for( tie(e, e_end)=boost::in_edges(node, wgraph->graph); e!=e_end; ++e ) {
    // get the WorkPiece id number, the output that it supplies, and the input that receives it
    const unsigned int id = wgraph->graph[boost::source(*e, wgraph->graph)]->piece->ID();
    const unsigned int inNum = wgraph->graph[*e]->inputDim;
    const unsigned int outNum = wgraph->graph[*e]->outputDim;

    // try to find the WorkPiece in the other upstream nodes
    auto it = inMap.find(id);

    if( it==inMap.end() ) { // if we have not yet needed this WorkPiece ...
      // ... add it to the list and store the input/output pair
      inMap[id] = std::vector<std::pair<unsigned int, unsigned int> >(1, std::pair<unsigned int, unsigned int>(inNum, outNum));
    } else { // we have needed this WorkPiece
      // ... add the input/output pair
      inMap[id].push_back(std::pair<unsigned int, unsigned int>(inNum, outNum));
    }
  }

  return inMap;
}

void ModGraphPiece::FillOutputMap() {
  // clear the map
  valMap.clear();

  // loop over the run order
  for( auto it : runOrder ) {
    // the inputs to this WorkPiece
    const ref_vector<Eigen::VectorXd>& ins = GetNodeInputs(it);

    // evaluate the current map and store a vector of references to its output
    auto wPiece = wgraph->GetPiece(it);
    auto currPiece = std::dynamic_pointer_cast<ModPiece>(wPiece);

    if(!currPiece){
      // If it can't be cast to a ModPiece, check to see if the output can be cast to an Eigen vector
      ref_vector<Eigen::VectorXd> output;

      std::vector<boost::any> anyIns(ins.size());
      for(int i=0; i<ins.size(); ++i)
        anyIns.at(i) = boost::any(ins.at(i));

      std::vector<boost::any> const& anyOut = wPiece->Evaluate(anyIns);


      for(int i=0; i<anyOut.size(); ++i){
        Eigen::VectorXd const& temp = AnyConstCast(anyOut.at(i));
        output.push_back(std::cref(temp));
      }
      valMap[wgraph->GetPiece(it)->ID()] = output;

    }else{
      assert(currPiece);
      valMap[wgraph->GetPiece(it)->ID()] = ToRefVector(currPiece->Evaluate(ins));
    }
  }
}

ref_vector<Eigen::VectorXd> ModGraphPiece::GetNodeInputs(boost::graph_traits<Graph>::vertex_descriptor node) const {

  // how many inputs does this node require?
  const int numIns = wgraph->GetPiece(node)->numInputs;

  // get the inputs for this node
  const std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > >& inMap = InputNodes(node);

  Eigen::VectorXd empty;
  ref_vector<Eigen::VectorXd> ins(numIns, std::cref(empty));

  // loop through the edges again, now we know which outputs supply which inputs
  for( auto edge : inMap ) {
    // loop over the input/output pairs supplied by this input
    for( auto in_out : edge.second ) {
      ins[in_out.first] = valMap.at(edge.first)[in_out.second];
    }
  }

  return ins;
}

std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > ModGraphPiece::RequiredOutputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn, unsigned int const wrtOut) const {
  // the ID of the current node
  const unsigned int nodeID = filtered_graphs[wrtIn]->operator[](node)->piece->ID();

  // get the outputs of this node that depend on the specified input
  std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > requiredOuts;

  if( nodeID==outputID ) { // if it is the output node ...
    // ... the user specifies the output derivative
    requiredOuts.push_back(std::tuple<unsigned int, unsigned int, unsigned int>(nodeID, wrtOut, wrtIn));

    return requiredOuts;
  }

  // loop though the output nodes
  boost::graph_traits<FilteredGraph>::out_edge_iterator eout, eout_end;
  for( tie(eout, eout_end)=boost::out_edges(node, *filtered_graphs[wrtIn]); eout!=eout_end; ++eout ) {
    // get the output number
    const unsigned int id = wgraph->GetPiece(boost::target(*eout, *filtered_graphs[wrtIn]))->ID();
    const unsigned int outNum = filtered_graphs[wrtIn]->operator[](*eout)->outputDim;
    const unsigned int inNum = filtered_graphs[wrtIn]->operator[](*eout)->inputDim;

    // if we have not already required this output, save it
    auto it = std::find(requiredOuts.begin(), requiredOuts.end(), std::tuple<unsigned int, unsigned int, unsigned int>(id, outNum, inNum));
    if( it==requiredOuts.end() ) {
      requiredOuts.push_back(std::tuple<unsigned int, unsigned int, unsigned int>(id, outNum, inNum));
    }
  }

  return requiredOuts;
}

std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > ModGraphPiece::RequiredInputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn) const {
  // how many inputs does this node require?
  const int numIns = filtered_graphs[wrtIn]->operator[](node)->piece->numInputs;

  std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > requiredIns;
  requiredIns.reserve(numIns);

  // loop though the output nodes
  boost::graph_traits<FilteredGraph>::in_edge_iterator ein, ein_end;
  for( tie(ein, ein_end)=boost::in_edges(node, *filtered_graphs[wrtIn]); ein!=ein_end; ++ein ) {
    // get the WorkPiece id number, the output that it supplies, and the input that receives it
    const unsigned int id = wgraph->GetPiece(boost::source(*ein, *filtered_graphs[wrtIn]))->ID();
    const unsigned int outNum = filtered_graphs[wrtIn]->operator[](*ein)->outputDim;
    const unsigned int inNum = filtered_graphs[wrtIn]->operator[](*ein)->inputDim;

    // store the requried input
    requiredIns.push_back(std::tuple<unsigned int, unsigned int, unsigned int>(id, outNum, inNum));
  }

  return requiredIns;
}

std::vector<int> ModGraphPiece::MatchInputs(std::shared_ptr<ModGraphPiece> otherPiece) const
{

  std::vector<std::shared_ptr<ConstantVector>> otherIns = otherPiece->GetConstantPieces();
  std::shared_ptr<WorkGraph> otherGraph = otherPiece->GetGraph();

  std::vector<int> outputs(otherIns.size());

  for(int i=0; i<otherIns.size(); ++i)
  {
    // get the downstream node and input index corresponding to this constant piece
    std::string constName = otherGraph->GetName( otherIns.at(i) );
    std::string sharedName = otherGraph->GetChildren( constName ).at(0);

    // Now try to find the same node and input in *this graph
    if(wgraph->HasNode(sharedName)){
      int inputIndex = otherGraph->GetEdges( constName, sharedName ).at(0).second;

      std::string upstreamName = wgraph->GetParent( sharedName, inputIndex);
      assert(upstreamName.size()>0);

      // Now, get the parent piece and check it against all of the constant pieces
      auto iter = std::find(constantPieces.begin(), constantPieces.end(), wgraph->GetPiece(upstreamName));
      if(iter != constantPieces.end()){
        outputs.at(i) = std::distance(constantPieces.begin(), iter);
      }else{
        outputs.at(i) = -1;
      }

    }else{
      outputs.at(i) = -1;
    }
  }

  return outputs;
}
