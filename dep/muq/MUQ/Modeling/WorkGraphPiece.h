#ifndef WORKGRAPHPIECE_H_
#define WORKGRAPHPIECE_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

#include "MUQ/Modeling/WorkPiece.h"
#include "MUQ/Modeling/ConstantPiece.h"
#include "MUQ/Modeling/NodeNameFinder.h"

namespace muq {
  namespace Modeling {

    /// This class keeps track of which nodes are downstream of a specified input
    class DependentPredicate {
    public:
      /// Required default constructor
      DependentPredicate();

      /**
	 @param[in] baseNode The input node --- we what to know which nodes are downstream of this node
	 @param[in] graph The graph holding all the nodes
      */
      DependentPredicate(boost::graph_traits<Graph>::vertex_descriptor const& baseNode, Graph const& graph);

      /**
	 @param[in] node Any node in the graph
	 \return true: the node is downstream of the graph, false: the node is not downstream of the graph
      */
      bool operator()(const boost::graph_traits<Graph>::vertex_descriptor& node) const;

    private:
      /// A vector of all the nodes downstream of the input node
      std::vector<boost::graph_traits<Graph>::vertex_descriptor> doesDepend;

      /**
	 @param[in] baseNode A node that depends on the input node (possible the input node itself)
	 @param[in] graph The graph holding all the nodes
      */
      void DownstreamNodes(const boost::graph_traits<Graph>::vertex_descriptor& baseNode, Graph const& graph);
    };

    /// Determine if the source of an edge is downstream of an input
    class DependentEdgePredicate {
    public:
      /// Required default constructor
      DependentEdgePredicate();

      /**
	 @param[in] nodePred All the downstream nodes of a given input
	 @param[in] graph The graph holding all of the nodes
      */
      DependentEdgePredicate(DependentPredicate nodePred, Graph const& graph);

      /**
	 @param[in] edge An edge in the graph
	 \return true: The source of the edge is a downstream node of the input, false: it is not a downstream node
      */
      bool operator()(const boost::graph_traits<Graph>::edge_descriptor& edge) const;

    private:
      /// The nodes that are downstream of the input
      DependentPredicate nodePred;

      /// The graph holding all the nodes
      const Graph* graph;
    };

    /// A filtered graph that only has nodes downstream of a specified input
    typedef boost::filtered_graph<Graph, DependentEdgePredicate, DependentPredicate> FilteredGraph;

    /// A muq::Modeling::WorkPiece created from a muq::Modeling::WorkGraph
    class WorkGraphPiece : public WorkPiece {
    public:

      /// Construct a muq::Modeling::WorkGraphPiece
      /**
	 Typically muq::Modeling::WorkGraphPiece's are constructed by calling muq::Modeling::WorkGraph::CreateWorkPiece
	 @param[in] graph From inputs to the output-of-interest
	 @param[in] constantPieces Pointers to the muq::Modeling::ConstantPiece's that hold the graph's inputs
	 @param[in] inputNames The names of each node in the graph corresponding to the constantPieces
	 @param[in] inTypes The input types (if known)
	 @param[in] outputNode The muq::Modeling::WorkPiece that we ultimately want to evaluate
	 @param[in] algebra Algebra to preform basic operations between different types (defaults to base class, which has common types)
       */
      WorkGraphPiece(std::shared_ptr<WorkGraph>                          wgraph,
                     std::vector<std::shared_ptr<ConstantPiece> > const& constantPieces,
                     std::vector<std::string>                     const& inputNames,
                     std::map<unsigned int, std::string>          const& inTypes,
                     std::shared_ptr<WorkPiece>                          outputNode);

      /// Default destructor
      virtual ~WorkGraphPiece();


      std::shared_ptr<WorkGraph> GetGraph(){return wgraph;};
      
    private:

      /// Evaluate each muq::Modeling::WorkPiece in the graph
      /**
	 @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
       */
      virtual void EvaluateImpl(ref_vector<boost::any> const& inputs) override;

   //    /// Compute the Jacobian for this muq::Modeling::WorkGraphPiece using the chain rule
   //    /**
	 // @param[in] wrtIn We are taking the Jacobian with respect to this input
	 // @param[in] wrtOut We are taking the Jacobian of this output
	 // @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
   //     */
   //    virtual void JacobianImpl(unsigned int const wrtIn, unsigned int const wrtOut, ref_vector<boost::any> const& inputs) override;
   //
   //    /// Compute the action of the Jacobian for this muq::Modeling::WorkGraphPiece using the chain rule
   //    /**
	 // @param[in] wrtIn We are taking the Jacobian with respect to this input
	 // @param[in] wrtOut We are taking the Jacobian of this output
	 // @param[in] vec We are applying the Jacobian to this object
	 // @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
   //     */
   //    virtual void JacobianActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override;
   //
   //    /// Compute the action of the Jacobian transpose for this muq::Modeling::WorkGraphPiece using the chain rule
   //    /**
	 // @param[in] wrtIn We are taking the Jacobian with respect to this input
	 // @param[in] wrtOut We are taking the Jacobian of this output
	 // @param[in] vec We are applying the Jacobian transpose to this object
	 // @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
   //     */
   //    virtual void JacobianTransposeActionImpl(unsigned int const wrtIn, unsigned int const wrtOut, boost::any const& vec, ref_vector<boost::any> const& inputs) override;

      /// Get the required outputs for a node in one of the filtered graphs
      /**
	 @param[in] node We want the outputs of this node
	 @param[in] wrtIn The input whose downstream nodes we care about
	 @param[in] wrtOut The output we are ultimately trying to differentiate wrt
	 \return The output nodes --- tuple: the output WorkPiece ID, the output number, and the input number
       */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > RequiredOutputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn, unsigned int wrtOut) const;

      /// Get the required inputs for a node in one of the filtered graphs
      /**
	 @param[in] node We want the inputs of this node
	 @param[in] wrtIn The input whose downstream nodes we care about
	 \return The input nodes --- tuple: the input WorkPiece ID, the output number, and the input number
       */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > RequiredInputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn) const;

      /// Fill the map from each node's muq::Modeling::WorkPiece::ID to its outputs
      void OutputMap();

      /// Set the inputs
      /**
	 Set the inputs in each the muq::Modeling::ConstantPiece.
       */
      void SetInputs(ref_vector<boost::any> const& inputs);

      /// Get the inputs from muq::Modeling::WorkGraphPiece::valMap to a specified node in the graph
      /**
	 @param[in] id The ID of the node of interest
	 \return A reference vector of inputs to that node
       */
      ref_vector<boost::any> Inputs(boost::graph_traits<Graph>::vertex_descriptor node) const;

      /// Get a the input nodes for a node
      /**
	 @param[in] node We want the input nodes for this node
	 \return A map from the input node's ID to the input/output number
      */
      std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > InputNodes(boost::graph_traits<Graph>::vertex_descriptor const& node) const;

      /// Run order computed during construction (input->output order)
      std::deque<boost::graph_traits<Graph>::vertex_descriptor> runOrder;

      // Run order for computing the derivatives of this muq::Modeling::WorkGraphPiece
      /**
	 Like muq::Modeling::WorkGraphPiece::runOrder, but specific to which input node is used (also in output->input order)
      */
      std::vector<std::deque<boost::graph_traits<Graph>::vertex_descriptor> > derivRunOrders;

      /// The WorkGraph associated with this WorkGraphPiece
      //std::shared_ptr<const Graph> graph;
      std::shared_ptr<WorkGraph> wgraph;

      std::vector<std::shared_ptr<FilteredGraph> > filtered_graphs;

      /// A the map from each node's muq::Modeling::WorkPiece::ID to its outputs
      std::unordered_map<unsigned int, ref_vector<boost::any> > valMap;

      /// The ID of the WorkPiece corresponding to the output node
      unsigned int outputID;

      /// The muq::Modeling::ConstantPiece's that store the inputs
      std::vector<std::shared_ptr<ConstantPiece> > constantPieces;

    };
  } // namespace Modeling
} // namespace muq

#endif
