#include "MUQ/Approximation/Polynomials/MonotoneExpansion.h"
#include "MUQ/Approximation/Polynomials/Monomial.h"

#include "MUQ/Approximation/Quadrature/GaussQuadrature.h"
#include "MUQ/Approximation/Polynomials/Legendre.h"

#include <memory>

using namespace muq::Approximation;
using namespace muq::Utilities;

MonotoneExpansion::MonotoneExpansion(std::shared_ptr<BasisExpansion> monotonePartsIn,
                                     bool                            coeffInput) : MonotoneExpansion({std::make_shared<BasisExpansion>(std::vector<std::shared_ptr<IndexedScalarBasis>>({std::make_shared<Monomial>()}))},
                                                                                                     {monotonePartsIn},
                                                                                                     coeffInput)
{
}

MonotoneExpansion::MonotoneExpansion(std::vector<std::shared_ptr<BasisExpansion>> const& generalPartsIn,
                                     std::vector<std::shared_ptr<BasisExpansion>> const& monotonePartsIn,
                                     bool                                                coeffInput) : ModPiece(GetInputSizes(monotonePartsIn, generalPartsIn, coeffInput),
                                                                                                                monotonePartsIn.size()*Eigen::VectorXi::Ones(1)),
                                                                                                       generalParts(generalPartsIn),
                                                                                                       monotoneParts(monotonePartsIn)
{
  assert(generalPartsIn.size() == monotonePartsIn.size());
  assert(generalPartsIn.at(0)->Multis()->GetMaxOrders().maxCoeff()==0);

  numInputs = -1;

  // Get the maximum order of of the monotone parts
  int maxOrder = 0;
  for(int i=0; i<monotoneParts.size(); ++i)
    maxOrder = std::max(maxOrder, monotoneParts.at(i)->Multis()->GetMaxOrders()(i) );

  // Number of quadrature points is based on the fact that an N point Gauss-quadrature
  // rule can integrate exactly a polynomial of order 2N-1
  int numQuadPts = ceil(0.5*(2.0*maxOrder + 1.0));
  GaussQuadrature gqSolver(std::make_shared<Legendre>(), numQuadPts-1);
  gqSolver.Compute(numQuadPts-1);

  quadPts = 0.5*(gqSolver.Points().transpose()+Eigen::VectorXd::Ones(numQuadPts));
  quadWeights = 0.5*gqSolver.Weights();
}

std::shared_ptr<MonotoneExpansion> MonotoneExpansion::Head(int numRows) const
{
  assert(numRows<=generalParts.size());

  std::vector<std::shared_ptr<BasisExpansion>> newGenerals(numRows);
  std::vector<std::shared_ptr<BasisExpansion>> newMonotones(numRows);

  for(int i=0; i<numRows; ++i){
    newGenerals.at(i) = generalParts.at(i);
    newMonotones.at(i) = monotoneParts.at(i);
  }

  return std::make_shared<MonotoneExpansion>(newGenerals, newMonotones);
}

Eigen::VectorXd MonotoneExpansion::EvaluateInverse(Eigen::VectorXd const& refPt) const{
  Eigen::VectorXd tgtPt0 = Eigen::VectorXd::Zero(refPt.size());
  return EvaluateInverse(refPt, tgtPt0);
}

Eigen::VectorXd MonotoneExpansion::EvaluateInverse(Eigen::VectorXd const& refPt,
                                                   Eigen::VectorXd const& tgtPt0) const{

  Eigen::VectorXd tgtPt = tgtPt0;
  Eigen::VectorXd mappedPt, resid, step, newPt;
  double residNorm = 1e4;
  double newResidNorm;
  const int maxLineIts = 20;

  while(residNorm>1e-11){
    resid = EvaluateForward(tgtPt) - refPt;
    step = JacobianWrtX(tgtPt).triangularView<Eigen::Lower>().solve(resid);

    double stepSize = 1.0;
    newPt = tgtPt - stepSize * step;

    for(int lineIt=0; lineIt<maxLineIts; ++lineIt){
      newResidNorm = (EvaluateForward(newPt) - refPt).norm();
      if(newResidNorm<(1.0-1e-9)*residNorm)
        break;
      stepSize *= 0.5;
      newPt = tgtPt - stepSize * step;
    }

    residNorm = newResidNorm;
    tgtPt = newPt;
  }

  return tgtPt;
}


unsigned MonotoneExpansion::NumTerms() const{

  unsigned numTerms = 0;
  for(int i=0; i<generalParts.size(); ++i)
    numTerms += generalParts.at(i)->NumTerms();
  for(int i=0; i<monotoneParts.size(); ++i)
    numTerms += monotoneParts.at(i)->NumTerms();

  return numTerms;
}

Eigen::VectorXd MonotoneExpansion::GetCoeffs() const
{
  Eigen::VectorXd output(NumTerms());

  unsigned currInd = 0;

  for(int i=0; i<generalParts.size(); ++i){
    output.segment(currInd, generalParts.at(i)->NumTerms()) = generalParts.at(i)->GetCoeffs().transpose();
    currInd += generalParts.at(i)->NumTerms();
  }

  for(int i=0; i<monotoneParts.size(); ++i){
    output.segment(currInd, monotoneParts.at(i)->NumTerms()) = monotoneParts.at(i)->GetCoeffs().transpose();
    currInd += monotoneParts.at(i)->NumTerms();
  }

  return output;
}

void MonotoneExpansion::SetCoeffs(Eigen::VectorXd const& allCoeffs)
{
  unsigned currInd = 0;

  for(int i=0; i<generalParts.size(); ++i){
    generalParts.at(i)->SetCoeffs(allCoeffs.segment(currInd, generalParts.at(i)->NumTerms()).transpose());
    currInd += generalParts.at(i)->NumTerms();
  }

  for(int i=0; i<monotoneParts.size(); ++i){
    monotoneParts.at(i)->SetCoeffs(allCoeffs.segment(currInd, monotoneParts.at(i)->NumTerms()).transpose());
    currInd += monotoneParts.at(i)->NumTerms();
  }
}

Eigen::VectorXd MonotoneExpansion::EvaluateForward(Eigen::VectorXd const& x) const{

  Eigen::VectorXd refPt = Eigen::VectorXd::Zero(x.size());

  // Fill in the general parts
  for(int i=0; i<generalParts.size(); ++i)
    refPt(i) += generalParts.at(i)->Evaluate(x).at(0)(0);

  // Now add the monotone part
  Eigen::VectorXd quadEvals(quadPts.size());
  for(int i=0; i<monotoneParts.size(); ++i){

    // Loop over quadrature points
    Eigen::VectorXd evalPt = x;//.head(i+1);
    evalPt(i) = x(i) * quadPts(0);

    double polyEval = monotoneParts.at(i)->Evaluate(evalPt).at(0)(0);

    quadEvals(0) = polyEval*polyEval;
    for(int k=1; k<quadPts.size(); ++k){
      evalPt(i) = x(i) * quadPts(k);
      polyEval = monotoneParts.at(i)->Evaluate(evalPt).at(0)(0);
      quadEvals(k) = polyEval*polyEval;
    }
    refPt(i) += quadEvals.dot(quadWeights)*x(i);
  }

  return refPt;
}
void MonotoneExpansion::EvaluateImpl(muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs)
{
  // Update the coefficients if need be
  if(inputs.size()>1){
    SetCoeffs(inputs.at(1).get());
  }

  outputs.resize(1);
  outputs.at(0) = EvaluateForward(inputs.at(0).get());
}


Eigen::MatrixXd MonotoneExpansion::JacobianWrtX(Eigen::VectorXd const& x) const{

  Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(x.size(), x.size());

  double polyEval;
  Eigen::MatrixXd polyGrad;

  // Add all of the general parts
  for(int i=0; i<generalParts.size(); ++i){
    Eigen::MatrixXd partialJac = generalParts.at(i)->Jacobian(0,0,x);
    jac.block(i,0,1,i) += partialJac.block(0,0,1,i);
  }

  // Add the monotone parts
  for(int i=0; i<monotoneParts.size(); ++i){
    Eigen::VectorXd evalPt = x;//.head(i+1);

    for(int k=0; k<quadPts.size(); ++k){
      evalPt(i) = x(i) * quadPts(k);
      polyEval = monotoneParts.at(i)->Evaluate(evalPt).at(0)(0);
      polyGrad = monotoneParts.at(i)->Jacobian(0,0,evalPt);
      jac.block(i,0,1,i) += 2.0 * x(i) * quadWeights(k) * polyEval * polyGrad.block(0,0,1,i);
      jac(i,i) += quadWeights(k) * (polyEval * polyEval + 2.0*x(i)*polyGrad(i)*polyEval*quadPts(k));
    }
  }

  return jac;
}


void MonotoneExpansion::JacobianImpl(unsigned int const                                wrtIn,
                                     unsigned int const                                wrtOut,
                                     muq::Modeling::ref_vector<Eigen::VectorXd> const& inputs)
{
  if(inputs.size()>1){
    SetCoeffs(inputs.at(1).get());
  }

  Eigen::VectorXd const& x = inputs.at(0);
  assert(x.size()==monotoneParts.size());

  Eigen::MatrixXd jac;
  if(wrtIn==0){

    jac = JacobianWrtX(x);

  }else if(wrtIn==1){

    // calculate the total number of coefficients
    unsigned numTerms = NumTerms();

    // Initialize the jacobian to zero
    jac = Eigen::MatrixXd::Zero(x.size(), numTerms);

    // Fill in the general portion of the jacobian
    unsigned currCoeff = 0;
    for(int i=0; i<generalParts.size(); ++i){
      Eigen::MatrixXd partialJac = generalParts.at(i)->Jacobian(0,1,x);
      jac.block(i, currCoeff, 1, generalParts.at(i)->NumTerms()) += partialJac;
      currCoeff += generalParts.at(i)->NumTerms();
    }

    // Fill in the monotone portion of the jacobian
    double polyEval;
    Eigen::MatrixXd polyGrad;
    for(int i=0; i<monotoneParts.size(); ++i){
      Eigen::VectorXd evalPt = x;//.head(i+1);

      for(int k=0; k<quadPts.size(); ++k){
        evalPt(i) = x(i) * quadPts(k);
        polyEval = monotoneParts.at(i)->Evaluate(evalPt).at(0)(0);
        polyGrad = monotoneParts.at(i)->Jacobian(0,1,evalPt);
        jac.block(i,currCoeff,1,monotoneParts.at(i)->NumTerms()) += 2.0 * x(i) * quadWeights(k) * polyEval * polyGrad;
      }

      currCoeff += monotoneParts.at(i)->NumTerms();
    }

  }

  jacobian = jac;
}


/** Returns the log determinant of the Jacobian (wrt x) matrix.
    Because the Jacobian is diagonal, the determinant is given by the
    product of the diagonal terms.  The log determinant is therefore
    the sum of the logs of the diagonal terms, which depend only on the
    monotone pieces of the parameterization.
*/
double MonotoneExpansion::LogDeterminant(Eigen::VectorXd const& evalPt,
                                         Eigen::VectorXd const& coeffs)
{
  SetCoeffs(coeffs);
  return LogDeterminant(evalPt);
}
double MonotoneExpansion::LogDeterminant(Eigen::VectorXd const& evalPt)
{
  assert(evalPt.rows() == inputSizes(0));
  muq::Modeling::ref_vector<Eigen::VectorXd> vecIn;
  vecIn.push_back(std::cref(evalPt));
  JacobianImpl(0,0,vecIn);
  return jacobian.diagonal().array().log().sum();
}

/** Returns the gradient of the Jacobian log determinant with respect to the coefficients. */
Eigen::VectorXd MonotoneExpansion::GradLogDeterminant(Eigen::VectorXd const& evalPt,
                                                      Eigen::VectorXd const& coeffs)
{
  SetCoeffs(coeffs);
  return GradLogDeterminant(evalPt);
}

Eigen::VectorXd MonotoneExpansion::GradLogDeterminant(Eigen::VectorXd const& x)
{

  Eigen::VectorXd output = Eigen::VectorXd::Zero(NumTerms());

  double polyEval;
  Eigen::MatrixXd polyGradX, polyGradC, polyD2;

  // Figure out what the first monotone coefficient is
  unsigned currInd = 0;
  for(int i=0; i<generalParts.size(); ++i)
    currInd += generalParts.at(i)->NumTerms();

  // Add the monotone parts to the gradient
  for(int i=0; i<monotoneParts.size(); ++i){
    Eigen::VectorXd evalPt = x;//.head(i+1);

    double part1 = 0.0;
    Eigen::VectorXd part2 = Eigen::VectorXd::Zero(monotoneParts.at(i)->NumTerms());

    for(int k=0; k<quadPts.size(); ++k){
      evalPt(i) = x(i) * quadPts(k);
      polyEval  = monotoneParts.at(i)->Evaluate(evalPt).at(0)(0);
      polyGradX = monotoneParts.at(i)->Jacobian(0,0,evalPt);
      polyGradC = monotoneParts.at(i)->Jacobian(0,1,evalPt);
      polyD2    = monotoneParts.at(i)->SecondDerivative(0, 0, 1, evalPt);

      part1 += quadWeights(k) * (polyEval * polyEval + 2.0*x(i)*polyGradX(i)*polyEval*quadPts(k));
      part2 += 2.0 * quadWeights(k) * (polyEval * polyGradC + x(i) * quadPts(k) * (polyEval*polyD2.row(i) + polyGradX(i)*polyGradC)).transpose();
    }

    output.segment(currInd, monotoneParts.at(i)->NumTerms()) = part2/part1;

    currInd += monotoneParts.at(i)->NumTerms();

  }

  return output;
}


Eigen::VectorXi MonotoneExpansion::GetInputSizes(std::vector<std::shared_ptr<BasisExpansion>> const& generalPartsIn,
                                                 std::vector<std::shared_ptr<BasisExpansion>> const& monotonePartsIn,
                                                 bool                                                coeffInput)
{
  Eigen::VectorXi output;
  if(coeffInput){

    int numCoeffs = 0;
    for(int i=0; i<generalPartsIn.size(); ++i)
      numCoeffs += generalPartsIn.at(i)->NumTerms();
    for(int i=0; i<monotonePartsIn.size(); ++i)
      numCoeffs += monotonePartsIn.at(i)->NumTerms();

    output.resize(2);
    output << monotonePartsIn.size(), numCoeffs;

  }else{
    output.resize(1);
    output << monotonePartsIn.size();
  }
  return output;
}
