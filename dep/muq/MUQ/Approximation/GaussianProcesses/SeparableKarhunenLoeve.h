#ifndef SEPARABLEKARHUNENLOEVE_H
#define SEPARABLEKARHUNENLOEVE_H

#include <Eigen/Core>

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"


namespace muq
{
namespace Approximation
{
    /**
    @class SeparableKarhunenLoeve
    @ingroup GaussianProcesses
    @seealso KarhunenLoeveExpansion
    @brief Implements KL expansions that take advantage of separable structure in both the domain and covariance kernel.
    @details

    ## Background ##
    Consider a two dimensional covariance kernel of the form \f$k((x,y),(x^\prime, y^\prime))\f$ that
    can be written as
    \f[
    k((x,y),(x^\prime, y^\prime)) = k_x(x,x^\prime)k_y(y,y^\prime).
    \f]
    This type of kernel is separable in all of its components.  It is also possible
    to have kernels that are partially separable.  For example, a three dimensional
    kernel \f$k_2((x,y,z),(x^\prime, y^\prime,z^\prime))\f$ might be separable
    in \f$x\f$ and \f$(y,z)\f$ but not all three components.  This would result in
    a decomposition of the form
    \f[
    k_2((x,y,z),(x^\prime, y^\prime,z^\prime)) = k_x(x,x^\prime) k_{yz}((y,z), (y^\prime,z^\prime)).
    \f]

    In both the separable and partially separable cases, it can often be more efficient
    to construct the KL decompositions for each component kernel (e.g., \f$k_x, k_y, \ldots\f$)
    seperately and then combine the individual expansions to form an expansion for the original
    kernel.  This is possible when the full-dimensional seed points are formed
    from a tensor product of points in the separable dimensions.
    */
    class SeparableKarhunenLoeve
    {

    public:

      /** Construct the separable expansion from a list of the seperate kernels
          and seed points/weights.

          @param[in] kernelsIn A vector of kernels for each separable component.
                               Note: Each individual kernel can be defined with multiple
                               input dimensions, but there can be no shared dimensions
                               between kernels.   As an example, the three component
                               kernel \f$k_2((x,y,z),(x^\prime, y^\prime,z^\prime))\f$
                               mentioned above would have one kernel that depended on
                               dimension "{0}" and another that depended on dimensions
                               "{1,2}", which is allowed.  However, the first kernel
                               could not depend on "{0,1}", because then both kernels
                               would depend on the second "y" component and the product
                               would no long be separable.
          @param[in] seedPtsIn A vector of seed points for each component.  The
                               length of this vector must match the length of the
                               kernels vector.  Each component of the vector is an \f$N\times M\f$
                               matrix, where \f$N\f$ is the number of inputs to the kernel (e.g., 1 for \f$k_x\f$ and 2 for \f$k_{yz}\f$)
                               and \f$M\f$ is the number of points used to define
                               the quadrature rule.
          @param[in] seedWtsIn A vector containing weights corresponding to the
                               seed points.  There must be the same number of wts
                               and pts for each component.
      */
      SeparableKarhunenLoeve(std::vector<std::shared_ptr<KernelBase>> kernelsIn,
                             std::vector<Eigen::MatrixXd> const& seedPtsIn,
                             std::vector<Eigen::VectorXd> const& seedWtsIn,
                             boost::property_tree::ptree options = boost::property_tree::ptree());

      virtual Eigen::MatrixXd GetModes(Eigen::Ref<const Eigen::MatrixXd> const& pts) const override;

      virtual unsigned int NumModes() const override;

    private:
      std::vector<std::shared_ptr<KarhunenLoeveBase>> components;

      std::shared_ptr<MultiIndexSet> modeInds;
      unsigned int numModes;

    }; // class SeparableKarhunenLoeve
  }
}


#endif
