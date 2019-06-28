project (UQTk)

SET(copy_FILES
  quadpnts.txt
  quadwghts.txt
  )
INSTALL(FILES ${copy_FILES} DESTINATION cpp/tests/QuadLUTest/)

add_executable (QuadLUTest main.cpp)

target_link_libraries (QuadLUTest uqtk  )

target_link_libraries (QuadLUTest depdsfmt )
target_link_libraries (QuadLUTest depcvode )
target_link_libraries (QuadLUTest depnvec  )
target_link_libraries (QuadLUTest depslatec)
target_link_libraries (QuadLUTest deplapack)
target_link_libraries (QuadLUTest depblas  )
target_link_libraries (QuadLUTest deplbfgs )
target_link_libraries (QuadLUTest depfigtree )
target_link_libraries (QuadLUTest depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (QuadLUTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (QuadLUTest ifcore)
  else()
    target_link_libraries (QuadLUTest ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (QuadLUTest gfortran stdc++)
  else()
    target_link_libraries (QuadLUTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

link_directories(${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/src/cvode)

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/mcmc  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/lapack)
include_directories(../../../dep/blas)
include_directories(../../../dep/slatec)
include_directories(../../../dep/lbfgs)
include_directories(../../../dep/cvode-2.7.0/include)

include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS QuadLUTest DESTINATION bin/tests/)

add_test(QuadLUTest QuadLUTest)
