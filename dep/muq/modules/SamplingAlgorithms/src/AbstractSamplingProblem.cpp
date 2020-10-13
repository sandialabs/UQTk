#include "MUQ/SamplingAlgorithms/AbstractSamplingProblem.h"

namespace muq{
  namespace SamplingAlgorithms{

    AbstractSamplingProblem::AbstractSamplingProblem
            (Eigen::VectorXi const& blockSizesIn, Eigen::VectorXi const& blockSizesQOIIn)
          : numBlocks(blockSizesIn.size()),
            blockSizes(blockSizesIn),
            numBlocksQOI(blockSizesQOIIn.size()),
            blockSizesQOI(blockSizesQOIIn)
            {assert(blockSizes.size()==numBlocks); assert(blockSizesQOI.size()==numBlocksQOI);}

    AbstractSamplingProblem::AbstractSamplingProblem(Eigen::VectorXi const& blockSizesIn) : AbstractSamplingProblem(blockSizesIn, Eigen::VectorXi::Zero(0)) {}

    std::shared_ptr<SamplingState> AbstractSamplingProblem::QOI() {
      return nullptr;
    }

  }
}
