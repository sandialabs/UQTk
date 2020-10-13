#include "MUQ/Modeling/Distributions/Gaussian.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

Gaussian::Gaussian(unsigned int dim,
                   InputMask    extraInputs) : GaussianBase(dim, GetExtraSizes(dim, extraInputs)),
                                               mode(ModeFromExtras(extraInputs)),
                                               inputTypes(extraInputs),
                                               covPrec(Eigen::VectorXd::Ones(dim))
{
  ComputeNormalization();
}

Gaussian::Gaussian(Eigen::VectorXd const& muIn,
                   InputMask              extraInputs) : GaussianBase(muIn, GetExtraSizes(muIn.size(), extraInputs)),
                                                         mode(ModeFromExtras(extraInputs)),
                                                         inputTypes(extraInputs),
                                                         covPrec(Eigen::VectorXd::Ones(muIn.size()))
{
  ComputeNormalization();
}


Gaussian::Gaussian(Eigen::VectorXd const& muIn,
                   Eigen::MatrixXd const& obj,
                   Gaussian::Mode         modeIn,
                   InputMask              extraInputs) : GaussianBase(muIn, GetExtraSizes(muIn.size(), extraInputs)),
                                                         mode(modeIn),
                                                         inputTypes(extraInputs),
                                                         covPrec(obj)
{
  CheckInputTypes(extraInputs, mode);
  assert(mean.rows()==covPrec.rows());
  if(covPrec.cols()>1)
    assert(mean.rows()==covPrec.cols());

  if(covPrec.cols()>1)
    sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();

  ComputeNormalization();
}

void Gaussian::CheckInputTypes(InputMask extraInputs, Mode mode)
{
  if( (((extraInputs & ExtraInputs::DiagCovariance)>0) || ((extraInputs & ExtraInputs::FullCovariance)>0)) && (mode == Mode::Precision))
    throw std::logic_error("Extra arguments passed to Gaussian constructor do not match the covariance mode.");
  if( (((extraInputs & ExtraInputs::DiagPrecision)>0) || ((extraInputs & ExtraInputs::FullPrecision)>0)) && (mode == Mode::Covariance))
    throw std::logic_error("Extra arguments passed to Gaussian constructor do not match the covariance mode.");
}

Gaussian::Mode Gaussian::ModeFromExtras(InputMask extraInputs)
{
  if( ((extraInputs & ExtraInputs::DiagCovariance)>0) || ((extraInputs & ExtraInputs::FullCovariance)>0)){
    return Gaussian::Mode::Covariance;
  }else{
    return Gaussian::Mode::Precision;
  }
}

Eigen::VectorXi Gaussian::GetExtraSizes(unsigned dim, InputMask extraInputs)
{
  Eigen::VectorXi output(2);
  int numExtras = 0;

  if((extraInputs & ExtraInputs::Mean)>0){
    output(0) = dim;
    numExtras++;
  }

  if( ((extraInputs & ExtraInputs::DiagCovariance)>0) || ((extraInputs & ExtraInputs::DiagPrecision)>0)){
    output(numExtras) = dim;
    numExtras++;
  }


  if( ((extraInputs & ExtraInputs::FullCovariance)>0) || ((extraInputs & ExtraInputs::FullPrecision)>0)){
    assert(numExtras<2);
    output(numExtras) = dim*dim;
    numExtras++;
  }

  return output.head(numExtras);
}


Eigen::MatrixXd Gaussian::GetCovariance() const
{
  if(mode==Mode::Covariance){
    if(covPrec.cols()==1){
      return covPrec.col(0).asDiagonal();
    }else{
      return covPrec;
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().matrix().asDiagonal();
    }else{
      return covPrec.selfadjointView<Eigen::Lower>().llt().solve(Eigen::MatrixXd::Identity(mean.rows(),mean.rows()));
    }
  }
}

Eigen::MatrixXd Gaussian::GetPrecision() const
{
  if(mode==Mode::Precision){
    if(covPrec.cols()==1){
      return covPrec.col(0).asDiagonal();
    }else{
      return covPrec;
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().matrix().asDiagonal();
    }else{
      return covPrec.selfadjointView<Eigen::Lower>().llt().solve(Eigen::MatrixXd::Identity(mean.rows(),mean.rows()));
    }
  }
}


void Gaussian::ComputeNormalization() {

  if( mode==Gaussian::Mode::Covariance ) {
    if(covPrec.cols()==1){
      logDet = covPrec.array().log().sum();
    }else{
      logDet = 2.0*std::log(sqrtCovPrec.matrixL().determinant());
    }
  } else if( mode==Gaussian::Mode::Precision ) {
    if(covPrec.cols()==1){
      logDet = -covPrec.array().log().sum();
    }else{
      logDet = -2.0*std::log(sqrtCovPrec.matrixL().determinant());
    }
  }
}

void Gaussian::ResetHyperparameters(ref_vector<Eigen::VectorXd> const& inputs)
{
  unsigned currInd = 0;

  if((inputTypes & ExtraInputs::Mean)>0){
    assert(inputs.at(currInd).get().size() == mean.size());
    mean = inputs.at(currInd);
    currInd++;
  }

  if(((inputTypes & ExtraInputs::DiagCovariance)>0) || ((inputTypes & ExtraInputs::DiagPrecision)>0)){

    if(inputs.at(currInd).get().size() != mean.size())
      throw muq::WrongSizeError("The given diagonal covariance or precision has " + std::to_string(inputs.at(currInd).get().size()) + " components, but " + std::to_string(mean.size()) + " were expected.");

    covPrec = inputs.at(currInd).get();
  }else if(((inputTypes & ExtraInputs::FullCovariance)>0) || ((inputTypes & ExtraInputs::FullPrecision)>0)){

    if(inputs.at(currInd).get().size() != mean.size()*mean.size())
      throw muq::WrongSizeError("The given covariance or precision has " + std::to_string(inputs.at(currInd).get().size()) + " components, but " + std::to_string(mean.size()*mean.size()) + " were expected.");

    Eigen::Map<const Eigen::MatrixXd> mat(inputs.at(currInd).get().data(), mean.size(), mean.size());
    covPrec = mat;
    sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();
  }

  ComputeNormalization();
}

void Gaussian::SetCovariance(Eigen::MatrixXd const& newCov) {

  mode = Gaussian::Mode::Covariance;

  assert(newCov.rows() == mean.rows());

  if(newCov.cols()>1)
    assert(newCov.cols() == mean.rows());

  covPrec = newCov;
  sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();

  // recompute the scaling constant
  ComputeNormalization();
}

void Gaussian::SetPrecision(Eigen::MatrixXd const& newPrec) {

  mode = Gaussian::Mode::Precision;

  assert(newPrec.rows() == mean.rows());

  if(newPrec.cols()>1)
    assert(newPrec.cols() == mean.rows());

  covPrec = newPrec;
  sqrtCovPrec = covPrec.selfadjointView<Eigen::Lower>().llt();

  // recompute the scaling constant
  ComputeNormalization();
}

Eigen::MatrixXd Gaussian::ApplyCovSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  if(mode==Gaussian::Mode::Covariance){
    if(covPrec.cols()==1){
      return covPrec.col(0).array().sqrt().matrix().asDiagonal()*x;
    }else{
      return sqrtCovPrec.matrixL()*x;
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().sqrt().matrix().asDiagonal() * x;
    }else{
      return sqrtCovPrec.matrixL().solve(x);
    }
  }
}
Eigen::MatrixXd Gaussian::ApplyPrecSqrt(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  if(mode==Gaussian::Mode::Precision){
    if(covPrec.cols()==1){
      return covPrec.col(0).array().sqrt().matrix().asDiagonal()*x;
    }else{
      return sqrtCovPrec.matrixL()*x;
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().sqrt().matrix().asDiagonal() * x;
    }else{
      return sqrtCovPrec.matrixL().solve(x);
    }
  }
}

Eigen::MatrixXd Gaussian::ApplyCovariance(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  if(mode==Gaussian::Mode::Covariance){
    if(covPrec.cols()==1){
      return covPrec.col(0).matrix().asDiagonal()*x;
    }else{
      return covPrec.selfadjointView<Eigen::Lower>()*x;
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().matrix().asDiagonal() * x;
    }else{
      return sqrtCovPrec.solve(x);
    }
  }
}
Eigen::MatrixXd Gaussian::ApplyPrecision(Eigen::Ref<const Eigen::MatrixXd> const& x) const
{
  if(mode==Gaussian::Mode::Covariance){
    if(covPrec.cols()==1){
      return covPrec.col(0).array().inverse().matrix().asDiagonal()*x;
    }else{
      return sqrtCovPrec.solve(x);
    }
  }else{
    if(covPrec.cols()==1){
      return covPrec.col(0).asDiagonal() * x;
    }else{
      return covPrec.selfadjointView<Eigen::Lower>()*x;
    }
  }
}


std::shared_ptr<Gaussian> Gaussian::Condition(Eigen::MatrixXd const& obsMat,
                                              Eigen::VectorXd const& data,
                                              Eigen::MatrixXd const& obsCov) const
{
  if(obsMat.cols() != Dimension())
    throw muq::WrongSizeError("In Gaussian::Condition, the number of columns in the observation matrix (" + std::to_string(obsMat.cols()) + ") does match the distribution dimension (" + std::to_string(Dimension()) + ")");
  if(obsMat.rows() != data.rows())
    throw muq::WrongSizeError("In Gaussian::Condition, the number of rows in the observation matrix (" + std::to_string(obsMat.rows()) + ") does match the size of the data (" + std::to_string(data.rows()) + ")");
  if(obsCov.rows() != data.rows())
    throw muq::WrongSizeError("In Gaussian::Condition, the length of the data vector (" + std::to_string(data.rows()) + ") does match the size of the observation covariance (" + std::to_string(obsCov.rows()) + ")");
  if((obsCov.rows() != obsCov.cols()) && (obsCov.cols()!=1))
    throw muq::WrongSizeError("In Gaussian::Condition, the given observation covariance has size " + std::to_string(obsCov.rows()) + "x" + std::to_string(obsCov.cols()) + " but should be square.");

  Eigen::MatrixXd HP = ApplyCovariance(obsMat.transpose()).transpose();
  Eigen::MatrixXd S = obsMat * HP.transpose();
  if(obsCov.cols()==1){
    S += obsCov.asDiagonal();
  }else{
    S += obsCov;
  }

  Eigen::MatrixXd K = S.selfadjointView<Eigen::Lower>().ldlt().solve(HP).transpose();

  Eigen::VectorXd postMu = mean + K*(data - obsMat*mean);
  Eigen::MatrixXd postCov = GetCovariance() - K*HP;

  return std::make_shared<Gaussian>(postMu, postCov);
}
