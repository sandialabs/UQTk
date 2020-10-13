#ifndef KERNELBASE_H
#define KERNELBASE_H


#include <assert.h>
#include <memory>
#include <iostream>
#include <fstream>

#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include "MUQ/Approximation/TemplatedArrayUtilities.h"

//#include "MUQ/Approximation/GaussianProcesses/StateSpaceGP.h"

#include "MUQ/Utilities/Exceptions.h"

#include <boost/property_tree/ptree.hpp>


namespace muq
{
namespace Modeling{
    class LinearSDE;
}
namespace Utilities{
    class LinearOperator;
}

namespace Approximation
{

/** @class KernelBase
    @ingroup CovarianceKernels
    @brief Base class for all covariance kernels.
*/
class KernelBase : public std::enable_shared_from_this<muq::Approximation::KernelBase>
{

public:

    KernelBase(unsigned int inputDimIn,
               unsigned int coDimIn,
               unsigned int numParamsIn) : KernelBase(inputDimIn, BuildDimInds(inputDimIn), coDimIn, numParamsIn)
    {};

    KernelBase(unsigned int              inputDimIn,
               std::vector<unsigned int> dimIndsIn,
               unsigned int              coDimIn,
               unsigned int              numParamsIn) : dimInds(dimIndsIn), inputDim(inputDimIn), coDim(coDimIn), numParams(numParamsIn)
    {
      assert(inputDim>0);
      assert(coDim>0);
    };


    virtual ~KernelBase(){};


    //virtual std::shared_ptr<muq::Approximation::KernelBase> GetPtr() {
    //    return shared_from_this();
    //}

    /// Overridden by ProductKernel
    virtual std::vector<std::shared_ptr<KernelBase>> GetSeperableComponents() {return std::vector<std::shared_ptr<KernelBase>>(1,Clone()); };


    virtual Eigen::MatrixXd Evaluate(Eigen::VectorXd const& x1,
                                     Eigen::VectorXd const& x2) const;

    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x) const;

    virtual Eigen::MatrixXd BuildCovariance(Eigen::MatrixXd const& x1,
                                            Eigen::MatrixXd const& x2) const;

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
                                Eigen::MatrixXd             const& ys,
                                Eigen::Ref<Eigen::MatrixXd>        cov) const;

    virtual void FillCovariance(Eigen::MatrixXd             const& xs,
                                Eigen::Ref<Eigen::MatrixXd>        cov) const;


    virtual void FillDerivCovariance(Eigen::MatrixXd             const& xs,
                                     Eigen::MatrixXd             const& ys,
                                     std::vector<int>            const& wrts,
                                     Eigen::Ref<Eigen::MatrixXd>        cov) const;

    /** @brief Returns derivatives of the kernel with respect to the first input, x1.
        @param[in] x1 The first position passed to the kernel.
        @param[in] x2 The second position passed to the kernel.
        @param[in] wrts A vector defining the order and directions of the
                   spatial derivatives.  wrts.size() is the derivative order.
        @return A matrix containing the derivatives of the kernel output with
                respect to the dimensions defined by wrts.
    */
    virtual Eigen::MatrixXd GetPosDerivative(Eigen::VectorXd  const& x1,
                                             Eigen::VectorXd  const& x2,
                                             std::vector<int> const& wrts) const;

    // virtual Eigen::MatrixXd GetParamDerivative(Eigen::VectorXd  const& x1,
    //                                            Eigen::VectorXd  const& x2,
    //                                            std::vector<int> const& dimWrts) const = 0;

    virtual Eigen::MatrixXd GetParamBounds() const
    {
      return paramBounds;
    };


    virtual Eigen::VectorXd GetParams() const{return cachedParams;};

    virtual void SetParams(Eigen::VectorXd const& params){assert(params.size()==numParams); cachedParams = params;};

    virtual std::shared_ptr<KernelBase> Clone() const = 0;


    /** @brief Returns a state space representation of the covariance kernel
        @details If this is a one dimensional kernel (i.e., inputDim=1 and coDim=1), this function returns a state space representation of the covariance kernel.  In particular, it returns a linear time invariant stochastic differential equation, whose solution, when started with the returned stationary covariance, provides the same information as this Gaussian process.   The first component of the vector-valued stochastic differential equation is related to the Gaussian process.  See "Kalman filtering and smoothing solutions to temporal Gaussian process regression models," by Jouni Hartikainen and Simo Sarkka, for more information.

    */
    virtual std::tuple<std::shared_ptr<muq::Modeling::LinearSDE>, std::shared_ptr<muq::Modeling::LinearOperator>, Eigen::MatrixXd> GetStateSpace(boost::property_tree::ptree sdeOptions=boost::property_tree::ptree()) const{
        throw muq::NotImplementedError("ERROR.  The GetStateSpace() function has not been implemented in this child of muq::Approximation::KernelBase.");
    };


    const std::vector<unsigned int> dimInds;


    const unsigned int inputDim;
    const unsigned int coDim;
    const unsigned int numParams;

    /** Evaluates a first or higher order derivative of the covariance kernel
        with respect to one of the position variables.
    */
    virtual void FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                                   Eigen::Ref<const Eigen::VectorXd> const& x2,
                                   Eigen::Ref<const Eigen::VectorXd> const& params,
                                   std::vector<int>                  const& wrts,
                                   Eigen::Ref<Eigen::MatrixXd>              block) const = 0;


    /** For particular points and parameters, this function fills in one
        block of the covariance matrix.
    */
    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const = 0;

protected:

    Eigen::VectorXd cachedParams;

    Eigen::MatrixXd paramBounds;

private:

    static std::vector<unsigned> BuildDimInds(unsigned dim)
    {
      std::vector<unsigned int> output(dim);
      for(unsigned int i=0; i<dim; ++i)
        output[i] = i;
      return output;
    }
};

}
}



#endif
