#ifndef MODGRAPHPIECE_H_
#define MODGRAPHPIECE_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

#include "MUQ/Modeling/ModPiece.h"
#include "MUQ/Modeling/WorkGraphPiece.h"
#include "MUQ/Modeling/ModGraphPiece.h"
#include "MUQ/Modeling/ConstantVector.h"
#include "MUQ/Modeling/NodeNameFinder.h"

namespace muq {
  namespace Modeling {

    /// A filtered graph that only has nodes downstream of a specified input
    typedef boost::filtered_graph<Graph, DependentEdgePredicate, DependentPredicate> FilteredGraph;

    /// A muq::Modeling::ModPiece created from a muq::Modeling::WorkGraph
    class ModGraphPiece : public ModPiece {
    public:

      /// Construct a muq::Modeling::ModGraphPiece
      /**
      	 Typically muq::Modeling::ModGraphPiece's are constructed by calling muq::Modeling::WorkGraph::CreateModPiece
      	 @param[in] graph From inputs to the output-of-interest
      	 @param[in] constantPieces Pointers to the muq::Modeling::ConstantVector's that hold the graph's inputs
      	 @param[in] inputNames The names of each node in the graph corresponding to the constantPieces
      	 @param[in] inTypes The input types (if known)
      	 @param[in] outputNode The muq::Modeling::WorkPiece that we ultimately want to evaluate
      	 @param[in] algebra Algebra to preform basic operations between different types (defaults to base class, which has common types)
       */
      ModGraphPiece(std::shared_ptr<WorkGraph>                           graph,
                    std::vector<std::shared_ptr<ConstantVector> > const& constantPieces,
                    std::vector<std::string>                      const& inputNames,
                    std::shared_ptr<ModPiece>                            outputNode);

      /// Default destructor
      virtual ~ModGraphPiece() = default;

      std::shared_ptr<WorkGraph> GetGraph(){return wgraph;};

      std::vector<std::shared_ptr<ConstantVector> > GetConstantPieces(){return constantPieces;};

      /**
        @brief Matches inputs with another ModGraphPiece by looking at node names and edges.

        @details Assume this ModGraphPiece has three inputs call \f$x_1\f$, \f$x_2\f$, and \f$y\f$,
        where \f$x_i\f$ denotes the \f$i^{th}\f$ input of a node \f$x\f$.  Let there
        be another ModGraphPiece with inputs \f$x_2\f$, \f$y\f$, and \f$z\f$.  This function
        tries to match the inputs between the two ModGraphPieces by looking at the node
        names and input indices.  It returns the input indices of *this that correspond to
        the inputs of the other piece, with the value "-1" reserved for nonoverlapping inputs.
        For the example above, the result would be \f$[1,2,-1]\f$.

        @params[in] otherPiece The other ModGraphPiece that we want to match with the inputs to this piece.
        @returns A std::vector containing indices into the inputs of this ModGraphPiece that correspond with the inputs of the otherPiece.  The length of this vector is equal to the number of inputs of otherPiece.
      */
      std::vector<int> MatchInputs(std::shared_ptr<ModGraphPiece> otherPiece) const;

    private:

      static Eigen::VectorXi ConstructInputSizes(std::vector<std::shared_ptr<ConstantVector> > const& constantPiecesIn);

      /// Evaluate each muq::Modeling::WorkPiece in the graph
      /**
	       @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
      */
      virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Compute the action of the Jacobian transpose for this muq::Modeling::WorkGraphPiece using the chain rule
      /**
      @param[in] wrtIn We are taking the Jacobian with respect to this input
      @param[in] wrtOut We are taking the Jacobian of this output
      @param[in] vec We are applying the Jacobian transpose to this object
      @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
       */
      // virtual void GradientImpl(unsigned int                const  outputDimWrt,
      //                           unsigned int                const  inputDimWrt,
      //                           ref_vector<Eigen::VectorXd> const& input,
      //                           Eigen::VectorXd             const& sensitivity);

      /// Compute the Jacobian for this muq::Modeling::WorkGraphPiece using the chain rule
      /**
         @param[in] wrtIn We are taking the Jacobian with respect to this input
         @param[in] wrtOut We are taking the Jacobian of this output
         @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
      */
      virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                                unsigned int                const  inputDimWrt,
                                ref_vector<Eigen::VectorXd> const& input) override;

      /// Compute the action of the Jacobian for this muq::Modeling::WorkGraphPiece using the chain rule
      /**
         @param[in] wrtIn We are taking the Jacobian with respect to this input
         @param[in] wrtOut We are taking the Jacobian of this output
         @param[in] vec We are applying the Jacobian to this object
         @param[in] inputs Inputs to the muq::Modeling::WorkGraphPiece
       */
      // virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
      //                                unsigned int                const  inputDimWrt,
      //                                ref_vector<Eigen::VectorXd> const& input,
      //                                Eigen::VectorXd             const& vec);

      /// Get the required outputs for a node in one of the filtered graphs
      /**
	       @param[in] node We want the outputs of this node
	       @param[in] wrtIn The input whose downstream nodes we care about
	       @param[in] wrtOut The output we are ultimately trying to differentiate wrt
	       @return The output nodes --- tuple: the output WorkPiece ID, the output number, and the input number
      */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > RequiredOutputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn, unsigned int wrtOut) const;

      /// Get the required inputs for a node in one of the filtered graphs
      /**
	       @param[in] node We want the inputs of this node
	       @param[in] wrtIn The input whose downstream nodes we care about
	       @return The input nodes --- tuple: the input WorkPiece ID, the output number, and the input number
       */
      std::vector<std::tuple<unsigned int, unsigned int, unsigned int> > RequiredInputs(boost::graph_traits<FilteredGraph>::vertex_descriptor const& node, unsigned int const wrtIn) const;

      /// Fill the map from each node's muq::Modeling::WorkPiece::ID to its outputs
      void FillOutputMap();

      /// Set the inputs
      /**
	       Set the inputs in each the muq::Modeling::ConstantPiece.
      */
      void SetInputs(ref_vector<Eigen::VectorXd> const& inputs);

      /// Get the inputs from muq::Modeling::WorkGraphPiece::valMap to a specified node in the graph
      /**
	       @param[in] id The ID of the node of interest
	       @return A reference vector of inputs to that node
       */
      ref_vector<Eigen::VectorXd> GetNodeInputs(boost::graph_traits<Graph>::vertex_descriptor node) const;

      /// Get a the input nodes for a node
      /**
	       @param[in] node We want the input nodes for this node
	       @return A map from the input node's ID to the input/output number
      */
      std::map<unsigned int, std::vector<std::pair<unsigned int, unsigned int> > > InputNodes(boost::graph_traits<Graph>::vertex_descriptor const& node) const;

      /// Run order computed during construction (input->output order)
      std::deque<boost::graph_traits<Graph>::vertex_descriptor> runOrder;

      // Run order for computing the derivatives of this muq::Modeling::WorkGraphPiece
      /**
	       Like muq::Modeling::WorkGraphPiece::runOrder, but specific to which input node is used (also in output->input order)
      */
      std::vector<std::deque<boost::graph_traits<Graph>::vertex_descriptor> > adjointRunOrders;

      /// The WorkGraph associated with this WorkGraphPiece
      std::shared_ptr<WorkGraph> wgraph;

      std::vector<std::shared_ptr<FilteredGraph> > filtered_graphs;

      /// A the map from each node's muq::Modeling::WorkPiece::ID to its outputs
      std::unordered_map<unsigned int, ref_vector<Eigen::VectorXd> > valMap;

      /// The ID of the WorkPiece corresponding to the output node
      unsigned int outputID;

      /// The muq::Modeling::ConstantVector's that store the inputs
      std::vector<std::shared_ptr<ConstantVector> > constantPieces;

    };
  } // namespace Modeling
} // namespace muq

#endif
