#include "MUQ/Approximation/Regression/Regression.h"

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#include "MUQ/Utilities/RandomGenerator.h"

#include "MUQ/Optimization/NLoptOptimizer.h"

#include "MUQ/Approximation/Polynomials/IndexedScalarBasis.h"

namespace pt = boost::property_tree;
using namespace muq::Utilities;
using namespace muq::Modeling;
using namespace muq::Optimization;
using namespace muq::Approximation;

Regression::Regression(pt::ptree const& pt) : WorkPiece(), order(pt.get<unsigned int>("Order")), inputDim(pt.get<unsigned int>("InputSize")), alpha(std::fmin(pt.get<double>("MaxPoisednessRadius", 1.0), 1.0)) {
  assert(alpha>0.0);
  poly = IndexedScalarBasis::Construct(pt.get<std::string>("PolynomialBasis", "Legendre"));

  // initalize the multi-index
  const std::string type = pt.get<std::string>("ExpansionType", "TotalOrder");
  if( type=="TotalOrder" ) {
    multi = MultiIndexFactory::CreateTotalOrder(inputDim, order);
  } else if( type=="Hyperbolic" ) {
    const double q = pt.get<double>("NormScale", 1.0);
    multi = MultiIndexFactory::CreateHyperbolic(inputDim, order, q);
  } else if( type=="Diagonal" ) {
    multi = std::make_shared<MultiIndexSet>(inputDim);

    // add a constant term
    std::shared_ptr<MultiIndex> constantIndex = std::make_shared<MultiIndex>(inputDim);
    multi->AddActive(constantIndex);

    // add the higher order terms
    for( unsigned int i=1; i<=order; ++i ) {
      for( unsigned int d=0; d<inputDim; ++d ) {
        auto singleIndex = MultiIndexFactory::CreateSingleTerm(inputDim, d, i);
        multi->AddActive(singleIndex);
      }
    }
  } else {
    std::cerr << std::endl;
    std::cerr << "ERROR: invalid polynomial expansion type in Regression.cpp" << std::endl;
    std::cerr << "\tChoose from: 'TotalOrder', 'Hyperbolic', or 'Diagonal'" << std::endl;
    std::cerr << std::endl;
    assert(false);
  }

  // set algorithm parameters for poisedness optimization with default values
  optPt.put("Ftol.AbsoluteTolerance", pt.get<double>("PoisednessConstant.Ftol.AbsoluteTolerance", 1.0e-8));
  optPt.put("Ftol.RelativeTolerance", pt.get<double>("PoisednessConstant.Ftol.RelativeTolerance", 1.0e-8));
  optPt.put("Xtol.AbsoluteTolerance", pt.get<double>("PoisednessConstant.Xtol.AbsoluteTolerance", 1.0e-8));
  optPt.put("Xtol.RelativeTolerance", pt.get<double>("PoisednessConstant.Xtol.RelativeTolerance", 1.0e-8));
  optPt.put("ConstraintTolerance", pt.get<double>("PoisednessConstant.ConstraintTolerance", 1.0e-8));
  optPt.put("MaxEvaluations", pt.get<unsigned int>("PoisednessConstant.MaxEvaluations", 1000));
  optPt.put("Algorithm", pt.get<std::string>("PoisednessConstant.Algorithm", "MMA"));
}

void Regression::EvaluateImpl(ref_vector<boost::any> const& inputs) {
  // if there are no points ... just return with an empty outputs
  if(inputs.size()==0) { return; }

  // make sure we can compute the Vandermonde matrix
  assert(multi);

  std::vector<Eigen::VectorXd> centered(inputs.size());
  for( unsigned int i=0; i<inputs.size(); ++i ) {
    // center and normalize
    centered[i] = (boost::any_cast<Eigen::VectorXd const>(inputs[i])-currentCenter).array()*currentRadius.inverse();
  }

  // get the Vandermonde matrix of the inputs
  const Eigen::MatrixXd vand = VandermondeMatrix(centered);
  assert(coeff.cols()==vand.cols());

  // compute the regression polynomial
  outputs.resize(1);
  outputs[0] = (Eigen::MatrixXd)(coeff*vand.transpose());
}

int Regression::NumInterpolationPoints() const {
  if( multi ) {
    return multi->Size();
  }

  std::cerr << std::endl << std::endl << "ERROR: Not able to compute the number of points required for interpolation" <<
    std::endl << "\tPolynomialRegressor.cpp NumInterpolationPoints()" << std::endl;
  assert(false);

  return -1;
}

void Regression::Fit(std::vector<Eigen::VectorXd> xs, std::vector<Eigen::VectorXd> const& ys, Eigen::VectorXd const& center) {
  assert(xs.size()>0);
  assert(xs.size()==ys.size());

  // set the current center
  currentCenter = center;

  // center the input points
  CenterPoints(xs);

  // compute basis coefficients
  coeff = ComputeCoefficients(xs, ys);
}

void Regression::Fit(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys) {
  assert(xs.size()>0);
  assert(xs.size()==ys.size());

  // preform the fit with zero center
  Fit(xs, ys, Eigen::VectorXd::Zero(xs[0].size()));
}

Eigen::MatrixXd Regression::ComputeCoefficients(std::vector<Eigen::VectorXd> const& xs, std::vector<Eigen::VectorXd> const& ys) const {
  assert(xs.size()==ys.size());

  // check to make sure we have more than the number of points required to interpolate
  const unsigned int interp = NumInterpolationPoints();
  if( xs.size()<interp ) {
    std::cerr << std::endl << "ERROR: Regression requires " << interp << " points to interpolate but only " << xs.size() << " are given." << std::endl;
    std::cerr << "\tTry fitting the regression with at least " << interp+1 << " points." << std::endl << std::endl;
    assert(xs.size()>NumInterpolationPoints());
  }

  // create the Vandermonde matrix and the rhs
  Eigen::MatrixXd vand = VandermondeMatrix(xs);
  const Eigen::MatrixXd rhs = ComputeCoefficientsRHS(vand, ys);
  vand = vand.transpose()*vand;

  // make the solver to do the regression
  auto solver = vand.colPivHouseholderQr();

  // comptue the coefficients
  return solver.solve(rhs).transpose();
}

Eigen::MatrixXd Regression::ComputeCoefficientsRHS(Eigen::MatrixXd const& vand, std::vector<Eigen::VectorXd> const& ys_data) const {
  // the dimension
  assert(ys_data.size()>0);
  const unsigned int dim = ys_data[0].size();

  // initialize space for the data
  Eigen::MatrixXd ys = Eigen::MatrixXd::Constant(ys_data.size(), dim, std::numeric_limits<double>::quiet_NaN());

  // copy the data into an Eigen type
  for( unsigned int i=0; i<ys_data.size(); ++i ) {
    ys.row(i) = ys_data[i];
  }

  // apply the Vandermonde matrix
  return vand.transpose()*ys;
}

Eigen::MatrixXd Regression::VandermondeMatrix(std::vector<Eigen::VectorXd> const& xs) const {
  assert(multi);

  // the number of points and the number of terms
  const unsigned int N = xs.size();
  const unsigned int M = multi->Size();
  assert(N>0);

  // initialize the matrix
  Eigen::MatrixXd vand = Eigen::MatrixXd::Ones(N, M);

  // each term is built by evaluating the polynomial basis
  for( unsigned int i=0; i<M; ++i ) { // loop through the terms
    // get the multi-index
    const Eigen::RowVectorXi alpha = multi->at(i)->GetVector();

    for( unsigned int pt=0; pt<N; ++pt ) { // loop through the points
      // get the point
      const Eigen::VectorXd& pnt = xs[pt];
      assert(alpha.size()==pnt.size());

      // each term is a product of 1D variables
      for( unsigned int v=0; v<alpha.size(); ++v ) {
	// evaluate the polynomial
	vand(pt, i) *= boost::any_cast<double const>(poly->Evaluate((unsigned int)alpha(v), pnt[v]) [0]);
      }
    }
  }

  return vand;
}

Eigen::ArrayXd Regression::CenterPoints(std::vector<Eigen::VectorXd>& xs, Eigen::VectorXd const& center) const {
  return CenterPoints(xs, center, xs.size());
}

Eigen::ArrayXd Regression::CenterPoints(std::vector<Eigen::VectorXd>& xs, Eigen::VectorXd const& center, unsigned int const kn) const {

  // is the center zero?
  const bool zeroCenter = center.norm()<std::numeric_limits<double>::epsilon();

  // loop through all of the input points and recenter them
  for( auto it=xs.begin(); it!=xs.end(); ++it ) {
    if( !zeroCenter ) {
      // recenter the the point
      *it -= center;
    }
  }

  // reset the radius
  Eigen::ArrayXd radius = Eigen::ArrayXd::Zero(center.size());
  for( auto it : xs ) { radius = radius.max(it.array().abs()); }

  for( auto it=xs.begin(); it!=xs.end(); ++it ) { (*it).array() *= radius.inverse(); }

  // return the radius
  return radius;
}

Eigen::ArrayXd Regression::CenterPoints(std::vector<Eigen::VectorXd>& xs) {
  // reset the current radius
  currentRadius = CenterPoints(xs, currentCenter, xs.size());
  return currentRadius;
}

std::pair<Eigen::VectorXd, double> Regression::PoisednessConstant(std::vector<Eigen::VectorXd> xs, Eigen::VectorXd const& center, int kn) const {
  // recenter so the points are on the unit ball (unit ball contains the first kn neighbors)
  assert(kn!=0);
  const unsigned int N = xs.size(); // the number of points
  if( kn<0 ) { kn = N; }
  assert(xs.size()>=kn);

  const Eigen::ArrayXd radius = CenterPoints(xs, center, kn);

  // the data for the lagrange polynomial are the Euclidean vectors (w.l.o.g. assume one dimensional output space)
  std::vector<Eigen::VectorXd> euclidean(N, Eigen::VectorXd::Zero(1));

  // a place to store the Lagrange coefficients
  std::vector<Eigen::RowVectorXd> lagrangeCoeff(N);

  // compute the coefficients for each Lagrange polynomial
  for( unsigned int i=0; i<N; ++i ) {
    if( i>0 ) { euclidean[i-1] = Eigen::VectorXd::Zero(1); }
    euclidean[i] = Eigen::VectorXd::Ones(1);

    // compute the coefficients
    lagrangeCoeff[i] = ComputeCoefficients(xs, euclidean);
  }

  auto cost = std::make_shared<PoisednessCost>(shared_from_this(), lagrangeCoeff, inputDim);
  auto constraint = std::make_shared<PoisednessConstraint>(inputDim, alpha);

  auto opt =
    std::make_shared<muq::Optimization::NLoptOptimizer>(cost, optPt);
  opt->AddInequalityConstraint(constraint);

  std::vector<Eigen::VectorXd> inputs;
  inputs.push_back((Eigen::VectorXd)Eigen::VectorXd::Zero(inputDim));

  const std::pair<Eigen::VectorXd, double>& soln =
    opt->Solve(inputs);
  assert(soln.second<=0.0); // we are minimizing a negative inner product so the optimal solution should be negative

  return std::pair<Eigen::VectorXd, double>((radius*soln.first.array()).matrix()+center, -soln.second);
}

void Regression::ComputeBasisDerivatives(Eigen::VectorXd const& point, std::vector<Eigen::VectorXd>& gradient) const {
  // the number of terms
  const unsigned int numCoeff = multi->Size();
  gradient.resize(numCoeff);

  // compute the gradient of the basis functions
  for( unsigned int i=0; i<numCoeff; ++i ) { // loop through the terms
    // initialize the gradient of basis function i
    gradient.at(i) = Eigen::VectorXd::Constant(inputDim, std::numeric_limits<double>::quiet_NaN());

    // the multiindex
    const Eigen::VectorXi multiIndex = multi->at(i)->GetVector();

    for( unsigned int d1=0; d1<inputDim; ++d1 ) { // loop through the dimensions of statespace
      // each term is a product of the 1D variables
      double result = 1.0;
      for( unsigned int v=0; v<inputDim; ++v ) {
	if( v==d1 ) { // the derivative we are differentiating w.r.t.
	  result *= poly->DerivativeEvaluate(multiIndex(v), 1, point(v));
	} else {
	  result *= boost::any_cast<double const>(poly->Evaluate((unsigned int)multiIndex(v), point(v)) [0]);
	}
      }

      // insert each entry into the vector
      gradient.at(i) (d1) = result;
    }
  }
}

Regression::PoisednessCost::PoisednessCost(std::shared_ptr<Regression const> parent, std::vector<Eigen::RowVectorXd> const& lagrangeCoeff, unsigned int const inDim) : CostFunction(Eigen::VectorXi::Constant(1, inDim)), parent(parent), lagrangeCoeff(lagrangeCoeff) {}

double Regression::PoisednessCost::CostImpl(ref_vector<Eigen::VectorXd> const& input) {
  const Eigen::VectorXd& x = input[0];

  const Eigen::RowVectorXd phi = parent->VandermondeMatrix(std::vector<Eigen::VectorXd>(1, x));

  Eigen::VectorXd lambda(lagrangeCoeff.size());
  for( unsigned int i=0; i<lagrangeCoeff.size(); ++i ) { lambda(i) = phi.dot(lagrangeCoeff[i]); }

  return -1.0*lambda.norm();
}

void Regression::PoisednessCost::GradientImpl(unsigned int const inputDimWrt, muq::Modeling::ref_vector<Eigen::VectorXd> const& input, Eigen::VectorXd const& sensitivity) {
  assert(inputDimWrt==0);
  const Eigen::VectorXd& x = input[0];

  const Eigen::RowVectorXd phi = parent->VandermondeMatrix(std::vector<Eigen::VectorXd>(1, x));

  // compute the gradient of the basis functions; each element is the gradient of a basis function
  std::vector<Eigen::VectorXd> gradBasis;
  parent->ComputeBasisDerivatives(x, gradBasis);
  assert(gradBasis.size()==lagrangeCoeff[0].size());
  assert(gradBasis[0].size()==inputSizes(0));

  gradient = Eigen::VectorXd::Zero(inputSizes(0));

  Eigen::VectorXd lambda(lagrangeCoeff.size());
  for( unsigned int i=0; i<lagrangeCoeff.size(); ++i ) {
    lambda(i) = phi.dot(lagrangeCoeff[i]);
    for( unsigned int j=0; j<gradBasis.size(); ++j ) {
      gradient += lambda(i)*lagrangeCoeff[i](j)*gradBasis[j];
    }
  }

  gradient *= -1.0*sensitivity(0)/lambda.norm();

  /*{ // sanity check ...
    std::cout << std::endl;
    Eigen::VectorXd gradFD = GradientByFD(0, 0, input, sensitivity);
    std::cout << "point (in GradientImpl): " << x.transpose() << std::endl;
    std::cout << "(cost) FD gradient: " << gradFD.transpose() << std::endl;
    std::cout << "(cost) gradient: " << gradient.transpose() << std::endl;
    }*/
}

Regression::PoisednessConstraint::PoisednessConstraint(unsigned int const inDim, double const alpha) :
  ModPiece(Eigen::VectorXi::Constant(1, inDim), Eigen::VectorXi::Ones(1)), alpha(alpha) {}

void Regression::PoisednessConstraint::EvaluateImpl(ref_vector<Eigen::VectorXd> const& input) {
  const Eigen::VectorXd& x = input[0];
  outputs.resize(outputSizes[0]);
  outputs[0] = (x.dot(x)-alpha)*Eigen::VectorXd::Ones(1);
}

void Regression::PoisednessConstraint::JacobianImpl(unsigned int const outputDimWrt,
                                                    unsigned int const inputDimWrt,
                                                    muq::Modeling::ref_vector<Eigen::VectorXd> const& input) {

  assert(inputDimWrt==0);
  const Eigen::VectorXd& x = input[0];

  jacobian.resize(inputSizes[0], outputSizes[0]);
  jacobian = 2.0*x;

  /*{ // sanity check ...
    std::cout << std::endl;
    Eigen::VectorXd gradFD = GradientByFD(0, 0, input, sensitivity);
    std::cout << "point (in GradientImpl): " << x.transpose() << std::endl;
    std::cout << "(constraint) FD gradient: " << gradFD.transpose() << std::endl;
    std::cout << "(constraint) gradient: " << gradient.transpose() << std::endl;
    }*/
}
