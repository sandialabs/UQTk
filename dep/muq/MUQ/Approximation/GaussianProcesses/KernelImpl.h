#ifndef KERNELIMPL_H
#define KERNELIMPL_H

#include "MUQ/Approximation/GaussianProcesses/KernelBase.h"

#include <stan/math/fwd/scal.hpp>

namespace muq
{
namespace Approximation
{

/**

\class KernelImpl
\ingroup CovarianceKernels
\brief Base class in CRTP pattern for covariance kernels
\details This class provides common functionality (such as computing Covariance matrices) for all covariance kernels.  It uses the curiously recurring template pattern and requires that child classes implement the following functions
- void EvaluateImpl(VectorType1, VectorType2, MatType)
- void GetDerivative(VectorType1, VectorType2, MatType)
- Eigen::VectorXd GetParams()
- void SetParams(Eigen::VectorXd)
*/
template<typename ChildType>
class KernelImpl : public KernelBase
{


public:

    KernelImpl(unsigned inputDimIn,
               unsigned coDimIn,
               unsigned numParamsIn) : KernelBase(inputDimIn, coDimIn, numParamsIn){};

    KernelImpl(unsigned              inputDimIn,
               std::vector<unsigned> dimIndsIn,
               unsigned              coDimIn,
               unsigned              numParamsIn) : KernelBase(inputDimIn, dimIndsIn, coDimIn, numParamsIn){};


    virtual ~KernelImpl(){};

    virtual std::shared_ptr<KernelBase> Clone() const override
    {
      return std::make_shared<ChildType>(static_cast<ChildType const &>(*this));
    }

    virtual void FillBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                           Eigen::Ref<const Eigen::VectorXd> const& x2,
                           Eigen::Ref<const Eigen::VectorXd> const& params,
                           Eigen::Ref<Eigen::MatrixXd>              block) const override
    {

      Eigen::VectorXd x1slice = GetSegment(x1);
      Eigen::VectorXd x2slice = GetSegment(x2);
      static_cast<const ChildType*>(this)->FillBlockImpl(Eigen::Ref<const Eigen::VectorXd>(x1slice),Eigen::Ref<const Eigen::VectorXd>(x2slice),params,block);
    };


    virtual void FillPosDerivBlock(Eigen::Ref<const Eigen::VectorXd> const& x1,
                                   Eigen::Ref<const Eigen::VectorXd> const& x2,
                                   Eigen::Ref<const Eigen::VectorXd> const& params,
                                   std::vector<int>                  const& wrts,
                                   Eigen::Ref<Eigen::MatrixXd>              block) const override
    {
      FillPosDerivBlockImpl(x1,x2,params,wrts,block);
    }

    void FillPosDerivBlockImpl(Eigen::Ref<const Eigen::VectorXd>  const& x1,
                               Eigen::Ref<const Eigen::VectorXd>  const& x2,
                               Eigen::Ref<const Eigen::VectorXd>  const& params,
                               std::vector<int>                   const& wrts,
                               Eigen::Ref<Eigen::MatrixXd>               block) const
    {
      // STAN will sometimes fail with third derivatives
      //assert(wrts.size()<3);

      if(wrts.size()==0){
        Eigen::VectorXd x1slice = GetSegment(x1);
        Eigen::VectorXd x2slice = GetSegment(x2);
        static_cast<const ChildType*>(this)->FillBlockImpl(Eigen::Ref<const Eigen::VectorXd>(x1slice), Eigen::Ref<const Eigen::VectorXd>(x2slice), params, block);
        return;
      }


      Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, 1> x1Temp(x1.size());
      Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, 1> x2Temp(x2.size());

      for(int i=0; i<x1.size(); ++i){
        x1Temp(i).val_ = x1(i);
        x1Temp(i).d_ = 0.0;
        x2Temp(i).val_ = x2(i);
        x2Temp(i).d_ = 0.0;
      }
      x1Temp(wrts.at(0)).d_ = 1.0;

      if(wrts.size()==1){
        Eigen::Matrix<stan::math::fvar<double>, Eigen::Dynamic, Eigen::Dynamic> cov(coDim, coDim);
        static_cast<const ChildType*>(this)->template FillBlockImpl<stan::math::fvar<double>,double, stan::math::fvar<double>>(x1Temp, x2Temp, params, cov);

        for(int j=0; j<coDim; ++j){
          for(int i=0; i<coDim; ++i){
            block(i,j) = cov(i,j).d_;
          }
        }
        return;
      }


      Eigen::Matrix<stan::math::fvar<stan::math::fvar<double>>, Eigen::Dynamic, 1> x1Temp2(x1.size());
      Eigen::Matrix<stan::math::fvar<stan::math::fvar<double>>, Eigen::Dynamic, 1> x2Temp2(x1.size());
      for(int i=0; i<x1.size(); ++i){
        x1Temp2(i).val_ = x1Temp(i);
        x1Temp2(i).d_ = 0.0;
        x2Temp2(i).val_ = x2Temp(i);
        x2Temp2(i).d_ = 0.0;
      }
      x1Temp2(wrts.at(1)).d_ = 1.0;

      if(wrts.size()==2){
        Eigen::Matrix<stan::math::fvar<stan::math::fvar<double>>, Eigen::Dynamic, Eigen::Dynamic> cov(coDim,coDim);
        static_cast<const ChildType*>(this)->template FillBlockImpl<stan::math::fvar<stan::math::fvar<double>>,double,stan::math::fvar<stan::math::fvar<double>>>(x1Temp2, x2Temp2, params, cov);

        for(int j=0; j<coDim; ++j){
          for(int i=0; i<coDim; ++i){
            block(i,j) = cov(i,j).d_.d_;
          }
        }

        return;
      }

      Eigen::Matrix<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>, Eigen::Dynamic, 1> x1Temp3(x1.size());
      Eigen::Matrix<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>, Eigen::Dynamic, 1> x2Temp3(x1.size());
      for(int i=0; i<x1.size(); ++i){
        x1Temp3(i).val_ = x1Temp2(i);
        x1Temp3(i).d_ = 0.0;
        x2Temp3(i).val_ = x2Temp2(i);
        x2Temp3(i).d_ = 0.0;
      }
      x1Temp3(wrts.at(2)).d_ = 1.0;

      if(wrts.size()==3){

        Eigen::Matrix<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>, Eigen::Dynamic, Eigen::Dynamic> cov(coDim,coDim);
        static_cast<const ChildType*>(this)->template FillBlockImpl<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>,double,stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>>(x1Temp3, x2Temp3, params, cov);

        for(int j=0; j<coDim; ++j){
          for(int i=0; i<coDim; ++i){
            block(i,j) = cov(i,j).d_.d_.d_;
          }
        }

        return;
      }

      Eigen::Matrix<stan::math::fvar<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>>, Eigen::Dynamic, 1> x1Temp4(x1.size());
      Eigen::Matrix<stan::math::fvar<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>>, Eigen::Dynamic, 1> x2Temp4(x1.size());
      for(int i=0; i<x1.size(); ++i){
        x1Temp4(i).val_ = x1Temp3(i);
        x1Temp4(i).d_ = 0.0;
        x2Temp4(i).val_ = x2Temp3(i);
        x2Temp4(i).d_ = 0.0;
      }
      x1Temp4(wrts.at(3)).d_ = 1.0;

      if(wrts.size()==4){

        Eigen::Matrix<stan::math::fvar<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>>, Eigen::Dynamic, Eigen::Dynamic> cov(coDim,coDim);
        static_cast<const ChildType*>(this)->template FillBlockImpl<stan::math::fvar<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>>,double,stan::math::fvar<stan::math::fvar<stan::math::fvar<stan::math::fvar<double>>>>>(x1Temp4, x2Temp4, params, cov);

        for(int j=0; j<coDim; ++j){
          for(int i=0; i<coDim; ++i){
            block(i,j) = cov(i,j).d_.d_.d_.d_;
          }
        }

        return;
      }

      assert(wrts.size()<=4);
    };

    template<typename ScalarType>
    Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> GetSegment(Eigen::Ref<const Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>> const& input) const
    {
      Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> output(dimInds.size());
      for(int i=0; i<dimInds.size(); ++i)
        output(i) = input(dimInds.at(i));

      return output;
    }
};


}
}

#endif
