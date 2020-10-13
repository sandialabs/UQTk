# make sure that Eigen supports the "Ref" command
set(CMAKE_REQUIRED_INCLUDES ${EIGEN3_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <Eigen/Core>
  void foo1(Eigen::Ref<Eigen::VectorXf> x){
  int temp = x.size();
  };
  int main(){
   Eigen::MatrixXf temp = Eigen::MatrixXf::Ones(3,3);
   foo1(temp.col(2));
   foo1(temp.col(0).head(2));
    return 0;
   }
  "
  EIGEN3_REF_COMPILES)

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <Eigen/Core>
    #include <Eigen/Eigenvalues>
    int main(){
     Eigen::VectorXd diag, subDiag;
     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
     solver.computeFromTridiagonal(diag,subDiag);
     return 0;
     }
    "
    EIGEN3_TRIDIAGONALEIGEN_COMPILES)

set(CMAKE_REQUIRED_INCLUDES ${EIGEN3_INCLUDE_DIR})
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <Eigen/Core>
	#include <unsupported/Eigen/FFT>
	#include <complex>

    int main(){
     Eigen::VectorXd temp1 = Eigen::VectorXd::Random(1024);
	 Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> temp2;
	 Eigen::FFT<double> fft;
	 fft.fwd(temp2,temp1);
      return 0;
     }
    "
    EIGEN3_FFT_COMPILES)

if(NOT EIGEN3_REF_COMPILES OR NOT EIGEN3_FFT_COMPILES OR NOT EIGEN3_TRIDIAGONALEIGEN_COMPILES)
	set(EIGEN3_TEST_FAIL 1)
	set(MUQ_USE_EIGEN3 OFF)
else()
	set(EIGEN3_TEST_FAIL 0)
	set(MUQ_USE_EIGEN3 ON)
endif()
