#ifndef SUMPIECE_H_
#define SUMPIECE_H_

#include "MUQ/Modeling/ModPiece.h"

namespace muq{
namespace Modeling{

  /** @class SumPiece
      @ingroup Modeling
      @brief Componentwise addition of two or more vectors.
  */
  class SumPiece : public ModPiece {
  public:


    /** @param[in] dim The dimension of the input vectors we want to add.
        @param[in] numInputs The number of input vectors.
    */
    SumPiece(unsigned int dim, unsigned int numInputs=2);

    virtual ~SumPiece() = default;

  protected:
    virtual void EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) override;

    virtual void GradientImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input,
                              Eigen::VectorXd             const& sensitivity) override;

    virtual void JacobianImpl(unsigned int                const  outputDimWrt,
                              unsigned int                const  inputDimWrt,
                              ref_vector<Eigen::VectorXd> const& input) override;

    virtual void ApplyJacobianImpl(unsigned int                const  outputDimWrt,
                                   unsigned int                const  inputDimWrt,
                                   ref_vector<Eigen::VectorXd> const& input,
                                   Eigen::VectorXd             const& vec) override;

    // virtual void ApplyHessianImpl(unsigned int                const  outWrt,
    //                               unsigned int                const  inWrt1,
    //                               unsigned int                const  inWrt2,
    //                               ref_vector<Eigen::VectorXd> const& input,
    //                               Eigen::VectorXd             const& sens,
    //                               Eigen::VectorXd             const& vec) override;
  }; // class SumPiece

}
}



#endif // #ifndef SumPiece_H_
