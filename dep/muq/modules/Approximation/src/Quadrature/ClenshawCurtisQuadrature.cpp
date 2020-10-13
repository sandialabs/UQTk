#include "MUQ/Approximation/Quadrature/ClenshawCurtisQuadrature.h"

using namespace muq::Approximation;


ClenshawCurtisQuadrature::ClenshawCurtisQuadrature(bool nestedIn) : Quadrature(1),
                                                                    nested(nestedIn) {}

unsigned int ClenshawCurtisQuadrature::Exactness(unsigned int index) const
{
  return IndexToNumPoints(index)-1;
}

unsigned int ClenshawCurtisQuadrature::IndexToNumPoints(unsigned int index) const
{
  if(nested){
    switch (index) {
      case 0:  return 1;
      case 1:  return 3;
      case 2:  return 5;
      case 3:  return 9;
      case 4:  return 17;
      case 5:  return 33;
      case 6:  return 65;
      case 7:  return 129;
      case 8:  return 257;
      case 9:  return 513;
      case 10: return 1025;
      case 11: return 2049;
      case 12: return 4097;
      case 13: return 8193;
      case 14: return 16385;
      case 15: return 32769;
      default: {
        std::stringstream msg;
        msg << "Requested a nested Clenshaw-Curtis rule with index " << index << ", which is not defined.  "
            << "The maximum nested index allowed is 15, which already has 32769 points.  "
            << "Do you really need more than that?" << std::endl;

        throw std::runtime_error(msg.str());
        return 0;
      }
    }
  }else{
    return index + 1;
  }
}

void ClenshawCurtisQuadrature::Compute(unsigned int index) {

  unsigned int numPoints  = IndexToNumPoints(index);

  pts.resize(1,numPoints);
  wts.resize(numPoints);

  if(numPoints==1){
    pts(0,0) = 0.0;
    wts(0) = 2.0;
    return;
  }

  for(int i=0; i<numPoints; ++i)
    pts(0,i) = std::cos( double(numPoints-i-1) * pi / double(numPoints - 1));

  pts(0,0) = -1.0;
  if((numPoints%2)==1)
    pts((numPoints+1)/2-1) = 0.0;
  pts(numPoints-1) = 1.0;

  wts = Eigen::VectorXd::Ones(numPoints);

  for(int i=0; i<numPoints; ++i) {
    double theta = double(i) * pi / double(numPoints-1);

    for(int j=0; j<(numPoints-1)/2; ++j){

      double b;
      if(2*(j+1)==(numPoints-1)){
        b = 1.0;
      }else{
        b = 2.0;
      }

      wts(i) -= b*std::cos(2.0*(j+1)*theta) / (4.0*std::pow(j+1,2) - 1);
    }
  }

  wts(0) /= double(numPoints-1);
  wts.segment(1,numPoints-2) *= 2.0/(numPoints-1);
  wts(numPoints-1) /= double(numPoints-1);

}
