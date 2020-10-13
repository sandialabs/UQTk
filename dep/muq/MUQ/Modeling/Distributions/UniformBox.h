#ifndef UNIFORMBOX_H_
#define UNIFORMBOX_H_

#include "MUQ/Utilities/VariadicMacros.h"
#include "MUQ/Utilities/Exceptions.h"

#include "MUQ/Modeling/Distributions/Distribution.h"


namespace muq {
  namespace Modeling {

    /** @class UniformBox
        @ingroup Distributions
        @brief Defines a normalized uniform distribution over a bounded rectangular domain.
@details This class defines a uniform distribution between coordinate-aligned lower and upper bounds.  Let \f$\theta \f$ denote the input random variable with components \f$\theta_i\f$ for \f$i\in \{1,2,\ldots, D\}\f$, let \f$L_i \f$ denote the lower bound for component \f$i\f$, and let \f$U_i\f$ denote the upper bound on \f$\theta_i \f$.  Then, this class implements \f[ \pi(\theta) = \prod \pi_i(\theta_i), \f] where \f[ \pi_i(\theta_i) = U\left[ L_i, U_i \right] = \left\{ \begin{array}{ll}\frac{1}{U_i-L_i} & L_i\leq x_i \leq U_i \\ 0 & \text{Otherwise} \end{array}\right. . \f]
This class can be initialized in three ways: with a matrix containing L_i and U_i, using std::pair<double,double> to collect the lower and upper bounds, and directly with a linear vector of doubles with entries that alternate between lower and upper bounds.  The following code snippets are all equivalent approaches for constructing the uniform distribution.

<h3>Construct from an Eigen::Matrix </h3>
@code{.cpp}
int dim = 4;

Eigen::MatrixXd bounds(dim,2);

bounds << 0.0, 1.0,
          0.5, 0.75,
          1.0, 2.0,
          4.0, 10.0;

auto dist = std::make_shared<UniformBox>(bounds);
@endcode

<h3>Construct from a std::vector of pairs </h3>
@code{.cpp}
int dim = 4;

std::vector<std::pair<double, double>> bounds(dim);

bounds.at(0) = std::make_pair(0.0, 1.0);
bounds.at(1) = std::make_pair(0.5, 0.75);
bounds.at(2) = std::make_pair(1.0, 2.0);
bounds.at(3) = std::make_pair(4.0, 10.0);

auto dist = std::make_shared<UniformBox>(bounds);
@endcode

<h3>Construct directly from lower and upper pairs</h3>
@code{.cpp}
int dim = 4;

std::pair<double, double> b1 = std::make_pair(0.0, 1.0);
std::pair<double, double> b2 = std::make_pair(0.5, 0.75);
std::pair<double, double> b3 = std::make_pair(1.0, 2.0);
std::pair<double, double> b4 = std::make_pair(4.0, 10.0);

auto dist = std::make_shared<UniformBox>(b1, b2, b3, b4);
@endcode

<h3> Construct with a linear vector of doubles </h3>
@code{.cpp}
int dim = 4;

std::vector<double> bounds(2*dim);

bounds.at(0) = 0.0;
bounds.at(1) = 1.0;
bounds.at(2) = 0.5;
bounds.at(3) = 0.75;
bounds.at(4) = 1.0;
bounds.at(5) = 2.0;
bounds.at(6) = 4.0;
bounds.at(7) = 10.0;

auto dist = std::make_shared<UniformBox>(bounds);
@endcode

<h3>Construct directly from doubles </h3>
@code{.cpp}
auto dist = std::make_shared<UniformBox>(0.0, 1.0, 0.5, 0.75, 1.0, 2.0, 4.0, 10.0);
@endcode
    */

    class UniformBox : public Distribution {
    public:

      virtual ~UniformBox() = default;

      /// Create a \f$d\f$-dimensional uniform distribution
      /**
	 @param[in] bounds A \f$d\f$-dimensional vector, each entry is a pair---first: lower bound in the dimension, second: upper bound in that dimension
       */
      UniformBox(Eigen::MatrixXd const& bounds);

      /// Create a \f$d\f$-dimensional uniform distribution
      /**
	 @param[in] args The upper/lower bound pairs (may be more than one)
       */
      template<typename... Args>
      inline UniformBox(std::pair<double,double> const& pair1, Args... args) : UniformBox(CreateBoundsPairs(pair1, args...)) {}

      inline UniformBox(std::vector<std::pair<double,double>> const& pairList) : UniformBox(CreateBoundsPairs(pairList)) {}

      template<typename... Args>
      inline UniformBox(double lb1, double ub1, Args... args) : UniformBox(CreateBoundsDouble(lb1, ub1, args...)) {}

      inline UniformBox(std::vector<double> const& bounds) : UniformBox(CreateBoundsDouble(bounds)) {}

    private:

      /// Evaluate the log-density
      /**
	 Inputs:
	 <ol>
	 <li> The state \f$x\f$
	 </ol>
	 \return The log-density (either 1 or negative infinity)
       */
      virtual double LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) override;

      /// Sample the distribution
      virtual Eigen::VectorXd SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) override;


      //static Eigen::MatrixXd CreateBounds(std::vector<double>& bounds);

      static Eigen::MatrixXd CreateBoundsPairs(std::vector<std::pair<double,double>> const& bounds);
      STATIC_VARIADIC_TO_VECTOR(CreateBoundsPairs, (std::pair<double, double>), (Eigen::MatrixXd))

      static Eigen::MatrixXd CreateBoundsDouble(std::vector<double> const& bounds);
      STATIC_VARIADIC_TO_VECTOR(CreateBoundsDouble, (double), (Eigen::MatrixXd))

          //static Eigen::MatrixXd CreateBounds(std::vector<std::pair<double,double>>& bounds);
          //STATIC_VARIADIC_TO_VECTOR_PART2(CreateBounds, double, (Eigen::MatrixXd))


      /// A matrix describing the bounding box.
      /**
	 A matrix that discribes the bounding box.  The first row contains the lower bounds and the second row corresponds to the upper bound.  Each column corresponds to a different inpu dimension.
       */
      const Eigen::MatrixXd bounds;

      // Use the values stored in bounds to compute the volume
      static double ComputeVolume(Eigen::MatrixXd const& boundsIn);

      // The volume of the region enclosed by the bounds
      const double volume;

    };
  } // namespace Modeling
} // namespace muq

#endif
