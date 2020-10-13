# make sure that the HDF5 library is available
set(CMAKE_REQUIRED_LIBRARIES ${SUNDIALS_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${SUNDIALS_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <idas/idas.h>

  #include <nvector/nvector_serial.h>

  #include <sundials/sundials_dense.h>
  #include <sundials/sundials_types.h>
  #include <sundials/sundials_math.h>

  static int idasJac(long int N,
                     realtype time,
  				   realtype alpha,
                     N_Vector state,
  				   N_Vector deriv,
                     N_Vector resid,
                     DlsMat   jac,
                     void    *user_data,
                     N_Vector tmp1,
                     N_Vector tmp2,
                     N_Vector tmp3){return 1;};

  int main()
  {
    N_Vector state;
    N_VDestroy(state);
    return 0;
  }


  "
  SUNDIALS_IDA_COMPILES)

if(NOT SUNDIALS_IDA_COMPILES)
  	set(SUNDIALS_TEST_FAIL 1)
else()
	  set(SUNDIALS_TEST_FAIL 0)
endif()
