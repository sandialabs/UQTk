#ifndef WORKGRAPHNODE_H_
#define WORKGRAPHNODE_H_

#include "MUQ/Modeling/WorkPiece.h"

namespace muq {
  namespace Modeling { 
    /// A node in a muq::Modeling::WorkGraph
    class WorkGraphNode {
    public:
      
      /// Create a node for muq::Modeling::WorkGraph
      /**
	 @param[in] piece A pointer to a muq::Modelling::Core::WorkPiece that will be called when this node is evaluated.
	 @input[in] name The name of the node
      */
      WorkGraphNode(std::shared_ptr<WorkPiece> piece, std::string const& name);
      
      /// A pointer to a muq::Modelling::WorkPiece that will be called when this node is evaluated.
      std::shared_ptr<WorkPiece> piece;
      
      /// The name of this node
      const std::string name;
      
    private:
      
    };
  } // namespace Modeling
} // namespace muq

#endif
