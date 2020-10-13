#include "MUQ/Approximation/PolynomialChaos/PolynomialChaosExpansion.h"

#include <fstream>
#include <iostream>

#include "MUQ/Approximation/Polynomials/OrthogonalPolynomial.h"

using namespace muq::Utilities;
using namespace muq::Approximation;

PolynomialChaosExpansion::PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet> const& multisIn,
                                                   Eigen::MatrixXd                                const& coeffsIn) : PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>>(multisIn->GetMultiLength(), basisCompsIn),
                                                                                                                                              multisIn,
                                                                                                                                              coeffsIn)
{}

PolynomialChaosExpansion::PolynomialChaosExpansion(std::shared_ptr<OrthogonalPolynomial>          const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet> const& multisIn,
                                                   unsigned int                                          outputDim) : PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>>(multisIn->GetMultiLength(), basisCompsIn),
                                                                                                                                               multisIn,
                                                                                                                                               outputDim)
{}

PolynomialChaosExpansion::PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet>   const& multisIn,
                                                   Eigen::MatrixXd                                  const& coeffsIn) : BasisExpansion(basisCompsIn, multisIn, coeffsIn)
{
}

PolynomialChaosExpansion::PolynomialChaosExpansion(std::vector<std::shared_ptr<IndexedScalarBasis>> const& basisCompsIn,
                                                   std::shared_ptr<muq::Utilities::MultiIndexSet>   const& multisIn,
                                                   unsigned int                                            outputDim) : BasisExpansion(basisCompsIn, multisIn, Eigen::MatrixXd::Zero(outputDim, multisIn->Size()))
{}



Eigen::VectorXd PolynomialChaosExpansion::GetNormalizationVec() const{

  Eigen::VectorXd result = Eigen::VectorXd::Zero(multis->Size());

  //compute each one
  for (unsigned int i = 0; i < multis->Size(); i++) {
    double norm = 1.0; //start the normalization at 1

    //loop over the dimensions and multiply the normalizations for each component polynomial
    Eigen::RowVectorXi multi = multis->IndexToMulti(i)->GetVector();
    for(int k=0; k<multi.size(); ++k){
      norm *= std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(k))->Normalization(multi(k));
    }
    result(i) = std::sqrt(norm);
  }

  return result;
};


Eigen::VectorXd PolynomialChaosExpansion::Variance() const
{
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //if there's only the one constant term, the PCE has variance zero
  if (normalVec.rows() <= 1) {
    return Eigen::VectorXd::Zero(coeffs.rows());
  }

  Eigen::VectorXd squareNorm = normalVec.tail(normalVec.rows() - 1).array().square();

  //for each output, the variance is the dot product of the squared coeffs with the squared norms of the PCE terms.
  //Thus, grab all but the leftmost column.
  //Have to normalize by what the constant integrates to.

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  //Since the variance is an expectation, we must normalize if the polynomials aren't quite set up
  //to integrate to one.
  return coeffs.rightCols(coeffs.cols() - 1).array().square().matrix() * squareNorm / dimMeasures.prod();
}

Eigen::MatrixXd PolynomialChaosExpansion::Covariance() const
{
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //if there's only the one constant term, the PCE has variance zero
  if (normalVec.rows() <= 1) {
    return Eigen::MatrixXd::Zero(coeffs.rows(),coeffs.rows());
  }

  Eigen::VectorXd squareNorm = normalVec.tail(normalVec.rows() - 1).array().square();

  //for each output, the variance is the dot product of the squared coeffs with the squared norms of the PCE terms.
  //Thus, grab all but the leftmost column.
  //Have to normalize by what the constant integrates to.

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  Eigen::VectorXd quadVec = squareNorm / dimMeasures.prod();

  // fill in upper portion of covariance matrix
  Eigen::MatrixXd cov(outputSizes(0), outputSizes(0));
  int npolys = coeffs.cols();
  for (int j = 0; j < outputSizes(0); ++j) {   // column index
    for (int i = j; i < outputSizes(0); ++i) { // row index
      cov(i,j) = coeffs.row(j).tail(npolys - 1) * quadVec.asDiagonal() * coeffs.row(i).tail(npolys - 1).transpose();
      cov(j,i) = cov(i,j);
    }
  }

  return cov;
}

Eigen::VectorXd PolynomialChaosExpansion::Mean() const
{
  return coeffs.col(0);
}


Eigen::VectorXd PolynomialChaosExpansion::Magnitude() const
{
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //take the matrix product between the squared coeffs and squared norms, then sqrt each term
  return (coeffs.array().square().matrix() * normalVec.array().square().matrix()).array().sqrt();
}



std::shared_ptr<PolynomialChaosExpansion> PolynomialChaosExpansion::ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>>    expansions,
                                                                                       Eigen::VectorXd                                    const& weights,
                                                                                       std::shared_ptr<MultiIndexSet>                     const& polynomials,
                                                                                       std::vector<std::vector<unsigned int>>             const& locToGlob)
{
  std::shared_ptr<PolynomialChaosExpansion> firstNonNull;
  for(unsigned int i=0; i<expansions.size();++i){
    if(expansions.at(i) != nullptr) {
      firstNonNull = expansions.at(i);
      break;
    }
  }
  assert(firstNonNull);

  // Loop over each expansion and add the weighted coefficients
  Eigen::MatrixXd newCoeffs = Eigen::MatrixXd::Zero(firstNonNull->outputSizes(0), polynomials->Size());

  for(unsigned int i=0; i<expansions.size(); ++i){
    if(std::abs(weights(i))>10.0*std::numeric_limits<double>::epsilon()){
      for(unsigned int k=0; k<expansions.at(i)->coeffs.cols(); ++k){
        newCoeffs.col(locToGlob.at(i).at(k)) += weights(i)*expansions.at(i)->coeffs.col(k);
      }
    }
  }

  return std::make_shared<PolynomialChaosExpansion>(firstNonNull->basisComps, polynomials, newCoeffs);
}

std::shared_ptr<PolynomialChaosExpansion> PolynomialChaosExpansion::ComputeWeightedSum(std::vector<std::shared_ptr<PolynomialChaosExpansion>> expansions,
                                                                                       Eigen::VectorXd                                 const& weights)
{
  assert(weights.size()==expansions.size());

  std::shared_ptr<PolynomialChaosExpansion> firstNonNull;
  for(unsigned int i=0; i<expansions.size();++i){
    if(expansions.at(i) != nullptr) {
      firstNonNull = expansions.at(i);
      break;
    }
  }
  assert(firstNonNull);

  unsigned int inputDim = firstNonNull->inputSizes(0);
  unsigned int outputDim = firstNonNull->outputSizes(0);

  // the collection of all polynomials that will be in the result
  auto allPolynomials = std::make_shared<MultiIndexSet>(inputDim);

  // Maps the local indices in each expansion to the global index of the sum
  std::vector<std::vector<unsigned int>> locToGlob(weights.size());

  // Make a union of all the multiindex sets, and keep track of the coefficient mapping
  for(int j=0; j<expansions.size(); ++j)
  {
    // if the weight is zero, skip this expansion
    if(std::abs(weights(j))>10.0*std::numeric_limits<double>::epsilon()){
      assert(expansions.at(j));

      //loop over all the polynomials in this expansion
      unsigned int numTerms = expansions.at(j)->multis->Size();
      locToGlob.at(j).resize(numTerms);
      for (unsigned int i = 0; i < numTerms; ++i)
        locToGlob.at(j).at(i) = allPolynomials->AddActive(expansions.at(j)->multis->IndexToMulti(i));
    }
  }

  return ComputeWeightedSum(expansions,weights,allPolynomials,locToGlob);
}


Eigen::VectorXd PolynomialChaosExpansion::TotalSensitivity(unsigned int targetDim) const
{

  //grab the total normalizations
  Eigen::VectorXd normalVec = GetNormalizationVec();

  // Zero out the terms that do not depend on the targetDim
  for (unsigned int i = 0; i < multis->Size(); ++i){
    if(multis->IndexToMulti(i)->GetValue(targetDim)==0)
      normalVec(i) = 0.0;
  }

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  //we want the sum of normalized involved coeffs divided by the normalized sum of all of them.
  //Don't forget the constant normalization the variance requires
  return (coeffs.array().square().matrix() * normalVec.array().square().matrix()).array() / dimMeasures.prod() / Variance().array();
}

Eigen::MatrixXd PolynomialChaosExpansion::TotalSensitivity() const
{
  Eigen::MatrixXd result(coeffs.rows(), inputSizes(0));

  for (unsigned int i = 0; i < inputSizes(0); ++i)
    result.col(i) = TotalSensitivity(i);

  return result;
}

Eigen::VectorXd PolynomialChaosExpansion::SobolSensitivity(unsigned int targetDim) const {
  return SobolSensitivity(std::vector<unsigned int>(1,targetDim));
}


Eigen::VectorXd PolynomialChaosExpansion::SobolSensitivity(std::vector<unsigned int> const& targetDims) const {
  //we're going to fill this index with either 1 or 0, depending on whether the term includes the target dimension
  Eigen::ArrayXd contributesToIndex = Eigen::VectorXd::Zero(NumTerms());

  // each term contributes to the main effects iff all nonzero terms in the multiindex correspond to a target dimension
  for(unsigned int i=0; i<multis->Size(); ++i){
    auto currMulti = multis->IndexToMulti(i);
    int miSum = currMulti->Sum();

    int partSum = 0;
    for(auto& k : targetDims)
      partSum += currMulti->GetValue(k);

    if((miSum>0)&&(partSum==miSum)){
      contributesToIndex(i) = 1;
    }
  }

  //grab the total normalizations
  Eigen::VectorXd normalVec = GetNormalizationVec();

  //filter the normalizations so that the ones that don't involve targetDim are zero
  Eigen::VectorXd filterdVec = normalVec.array() * contributesToIndex;

  Eigen::VectorXd dimMeasures(inputSizes(0));
  for (unsigned int i = 0; i < inputSizes(0); ++i)
    dimMeasures(i) = std::dynamic_pointer_cast<OrthogonalPolynomial>(basisComps.at(i))->Normalization(0);

  //we want the sum of normalized involved coeffs divided by the normalized sum of all of them.
  //Don't forget the constant normalization the variance requires
  return (coeffs.array().square().matrix() * filterdVec.array().square().matrix()).array() / dimMeasures.prod() / Variance().array();
}

Eigen::MatrixXd PolynomialChaosExpansion::MainSensitivity() const
{
  Eigen::MatrixXd result(coeffs.rows(), inputSizes(0));

  for (unsigned int i = 0; i < inputSizes(0); ++i)
    result.col(i) = SobolSensitivity(i);

  return result;
}
