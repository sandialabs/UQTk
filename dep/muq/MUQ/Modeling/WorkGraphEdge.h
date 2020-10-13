#ifndef WORKGRAPHEDGE_H_
#define WORKGRAPHEDGE_H_

namespace muq {
  namespace Modeling {
    /// An edge in a muq::Modeling::WorkGraph
    class WorkGraphEdge {
    public:

      /// Create an edge for muq::Modeling::WorkGraph
      /**
	      @param[in] outputDim The output dimension of the output node that will be given to the input node
	      @param[in] inputDim The input dimension of the input node that will be given the output of the output node
      */
      WorkGraphEdge(unsigned int const outputDim, unsigned int const inputDim);

      /// The output dimension
      const unsigned int outputDim;

      /// The input dimension
      const unsigned int inputDim;

    private:

    };
  } // namespace Modeling
} // namespace muq

#endif
