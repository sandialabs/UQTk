#ifndef DENSITYPRODUCT_H
#define DENSITYPRODUCT_H

#include "MUQ/Modeling/Distributions/Density.h"

namespace muq {
  namespace Modeling {

    class DensityProduct : public DensityBase {
    public:

      /** Constructs a density product accepting numPiecesIn inputs.  For example,
          a Bayesian posterior density is composed of a likelihood and prior, so
          numPiecesIn=2.
          @param[in] numPiecesIn Number of input densities.
      */
      DensityProduct(int numPiecesIn);

      virtual ~DensityProduct() = default;

    private:

      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      virtual Eigen::VectorXd GradLogDensityImpl(unsigned int wrt,
                                                 ref_vector<Eigen::VectorXd> const& inputs) override;

    };
  } // namespace Modeling
} // namespace muq

#endif
