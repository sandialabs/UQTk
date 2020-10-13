#ifndef EIGENUTILITIES_H
#define EIGENUTILITIES_H

#include <Eigen/Core>

namespace muq{
namespace Utilities{

  template<typename ScalarType>
  struct VectorLessThan
  {
      bool operator()(Eigen::Matrix<ScalarType, Eigen::Dynamic,1> const& left,
                      Eigen::Matrix<ScalarType, Eigen::Dynamic,1> const& right) const
      {
          if(left.size()!=right.size())
            return left.size()<right.size();

          for(unsigned int i=0; i<left.size(); ++i){
            if(left[i]<right[i]-5.0*std::numeric_limits<double>::epsilon()) return true;
            if(left[i]>right[i]+5.0*std::numeric_limits<double>::epsilon()) return false;
          }

          // At this point, they must be equal and therefore not less than
          return false;
      }
  };
}
}



#endif // #ifndef EIGENUTILITIES_H
