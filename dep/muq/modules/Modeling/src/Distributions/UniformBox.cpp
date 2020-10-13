#include "MUQ/Modeling/Distributions/UniformBox.h"

#include "MUQ/Utilities/RandomGenerator.h"

using namespace muq::Utilities;
using namespace muq::Modeling;

UniformBox::UniformBox(Eigen::MatrixXd const& boundsIn) : Distribution(boundsIn.rows(),
                                                                       Eigen::VectorXi()),
                                                          bounds(boundsIn),
                                                          volume(ComputeVolume(boundsIn)) {

}

double UniformBox::ComputeVolume(Eigen::MatrixXd const& boundsIn)
{
    return (boundsIn.col(1)-boundsIn.col(0)).prod();
}

Eigen::MatrixXd UniformBox::CreateBoundsPairs(std::vector<std::pair<double,double>> const& boundsVec)
{
    Eigen::MatrixXd bounds(boundsVec.size(),2);
    for(int i=0; i<boundsVec.size(); ++i){
        bounds(i,0) = boundsVec.at(i).first;
        bounds(i,1) = boundsVec.at(i).second;
    }

    return bounds;
}

Eigen::MatrixXd UniformBox::CreateBoundsDouble(std::vector<double> const& boundsVec)
{
    assert(boundsVec.size()%2==0);

    Eigen::MatrixXd bounds(boundsVec.size()/2,2);
    for(int i=0; i<boundsVec.size()/2; ++i){
        bounds(i,0) = boundsVec.at(2*i);
        bounds(i,1) = boundsVec.at(2*i+1);
    }

    return bounds;
}

double UniformBox::LogDensityImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  // loop though each dimension
  for( unsigned int i=0; i<bounds.rows(); ++i ) {
    // get the value of the point for this dimension
    const double x = inputs.at(0).get()(i);

    // if this dimension is outside the hypercube
    if( x<bounds(i,0) || x>bounds(i,1) ) {
      // return negative infinity
      return -std::numeric_limits<double>::infinity();
    }
  }

  // the point is inside the domain...
  return -log(volume);
}

Eigen::VectorXd UniformBox::SampleImpl(ref_vector<Eigen::VectorXd> const& inputs) {
  assert(bounds.rows()>0);

  Eigen::VectorXd sample;
  sample = bounds.col(0).array() + (bounds.col(1)-bounds.col(0)).array() * RandomGenerator::GetUniform(bounds.rows()).array();
  return sample;
}
