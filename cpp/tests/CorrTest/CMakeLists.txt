project (UQTk)

SET(copy_FILES
  set_1_w1.dat
  set_1_w2.dat
  set_2_w1.dat
  set_2_w2.dat
  set_3_w1.dat
  set_3_w2.dat
  set_4_w1.dat
  set_4_w2.dat
  set_5_w1.dat
  set_5_w2.dat
  set_6_w1.dat
  set_6_w2.dat
  )
file(COPY ${copy_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable (CorrTest main.cpp)

target_link_libraries (CorrTest uqtk )

target_link_libraries (CorrTest depdsfmt )
target_link_libraries (CorrTest depcvode )
target_link_libraries (CorrTest depnvec  )
target_link_libraries (CorrTest depslatec)
target_link_libraries (CorrTest deplapack)
target_link_libraries (CorrTest depblas  )
target_link_libraries (CorrTest deplbfgs  )
target_link_libraries (CorrTest depfigtree  )
target_link_libraries (CorrTest depann  )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (CorrTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (CorrTest ifcore ifport)
  else()
    target_link_libraries (CorrTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (CorrTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (CorrTest gfortran stdc++)
  else()
    target_link_libraries (CorrTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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
include_directories(../../../dep/figtree)
include_directories(../../../dep/cvode-2.7.0/include)

include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS CorrTest DESTINATION bin/tests/)

add_test(CorrTest CorrTest)
