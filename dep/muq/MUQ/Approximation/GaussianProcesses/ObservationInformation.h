#ifndef OBSERVATIONINFORMATION_H
#define OBSERVATIONINFORMATION_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"
#include "MUQ/Modeling/LinearAlgebra/LinearOperator.h"

#include <Eigen/Core>

#include <vector>

namespace muq{
namespace Approximation{

    /** @defgroup GP_Observations
        @ingroup GaussianProcesses
        @brief Tools for defining linear observations of Gaussian Processes.
    */

    /** @ingroup GP_Observations
        @class ObservationInformation
        @brief Class for defining linear observations of a Gaussian process
    */
    class ObservationInformation : public std::enable_shared_from_this<ObservationInformation>
    {
    public:

        ObservationInformation(std::shared_ptr<muq::Modeling::LinearOperator> Hin,
                               Eigen::Ref<const Eigen::VectorXd> const&       locIn,
                               Eigen::Ref<const Eigen::VectorXd> const&       obsIn,
                               Eigen::Ref<const Eigen::MatrixXd> const&       obsCovIn) : H(Hin), loc(locIn), obs(obsIn), obsCov(obsCovIn){};

        virtual ~ObservationInformation() = default;

        virtual void FillSelfCov(std::shared_ptr<KernelBase> kernel,
                                 Eigen::Ref<Eigen::MatrixXd> covBlock);

        virtual void FillCrossCov(Eigen::Ref<const Eigen::VectorXd> const& otherLoc,
                                  std::shared_ptr<KernelBase>              kernel,
                                  Eigen::Ref<Eigen::MatrixXd>              covBlock);

        virtual void FillCrossCov(std::shared_ptr<ObservationInformation> otherObs,
                                  std::shared_ptr<KernelBase>             kernel,
                                  Eigen::Ref<Eigen::MatrixXd>             covBlock);

        // The observation operator
        std::shared_ptr<muq::Modeling::LinearOperator> H;

        // The location of the observation
        Eigen::VectorXd loc;

        // The observed data
        Eigen::VectorXd obs;

        // The covariance of the observational noise
        Eigen::MatrixXd obsCov;

      protected:
        virtual Eigen::MatrixXd BuildBaseCovariance(Eigen::Ref<const Eigen::VectorXd> const& otherObs,
                                                    std::shared_ptr<KernelBase>              kernel);

        virtual Eigen::MatrixXd BuildBaseCovariance(std::shared_ptr<KernelBase>              kernel);

        virtual Eigen::MatrixXd BuildBaseCovariance(std::shared_ptr<ObservationInformation> otherObs,
                                                    std::shared_ptr<KernelBase>              kernel);

    };


    /** @ingroup GP_Observations
        @class DerivativeObservation
        @brief Class that defines an observation involving linear combinations of GP derivatives
        @details Let \f$y\f$ denote the observable random variable and let \f$u(x)\f$
        denote the Gaussian process.  This class defines observations of the form
        \f[
        y = H \left[ \begin{array}{c} \frac{ \partial^{N_1}u(x)}{\partial x_{n(1,1)} \ldots \partial x_{n(1,N_1)} }\\ \frac{ \partial^{N_2}u(x)}{\partial x_{n(2,1)} \ldots \partial x_{n(2,N_2)} } \\ \vdots \\ \frac{ \partial^{N_M}u(x)}{\partial x_{n(M,1)} \ldots \partial x_{n(M,N_M)} }  \end{array}\right],
        \f]
        for some appropriately sized matrix \f$H\f$.
    */
    class DerivativeObservation : public ObservationInformation
    {
    public:
        friend class ObservationInformation;

        DerivativeObservation(std::shared_ptr<muq::Modeling::LinearOperator> Hin,
                              Eigen::Ref<const Eigen::VectorXd> const&        locIn,
                              Eigen::Ref<const Eigen::VectorXd> const&        obsIn,
                              Eigen::Ref<const Eigen::MatrixXd> const&        obsCovIn,
                              std::vector<std::vector<int>>                   derivCoordsIn) : ObservationInformation(Hin, locIn, obsIn, obsCovIn),
                                                                                               derivCoords(derivCoordsIn){};

        virtual ~DerivativeObservation() = default;

        /**
        Derivatives to consider.  These define the \f$n(i,j)\f$ quantities
        in the description above.
        */
        std::vector<std::vector<int>> derivCoords;

      protected:

        virtual Eigen::MatrixXd BuildBaseCovariance(Eigen::Ref<const Eigen::VectorXd> const& otherObs,
                                                    std::shared_ptr<KernelBase>              kernel) override;

        virtual Eigen::MatrixXd BuildBaseCovariance(std::shared_ptr<KernelBase>              kernel) override;

        virtual Eigen::MatrixXd BuildBaseCovariance(std::shared_ptr<ObservationInformation> otherObs,
                                                    std::shared_ptr<KernelBase>              kernel) override;



    };


} // namespace muq
} // namespace Approximation





#endif // #ifndef OBSERVATIONINFORMATION_H
