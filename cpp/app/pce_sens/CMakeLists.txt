project (UQTk)

link_directories(${CMAKE_SUNDIALS_DIR}/lib)
link_directories (${INSTALL_LIB_DIR})
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/cvode)
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/nvector/serial)
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/sunlinsol/dense)
link_directories(/Users/kjohnst/UQTk-build/dep/sundials/src/sunmatrix/dense)

add_executable (pce_sens pce_sens.cpp)

target_link_libraries (pce_sens uqtkpce  )
target_link_libraries (pce_sens uqtkquad )
target_link_libraries (pce_sens uqtktools)
target_link_libraries (pce_sens uqtkarray)

target_link_libraries (pce_sens depdsfmt )
target_link_libraries (pce_sens sundials_cvode)
target_link_libraries (pce_sens sundials_nvecserial)
target_link_libraries (pce_sens sundials_sunlinsoldense)
target_link_libraries (pce_sens sundials_sunmatrixdense)
target_link_libraries (pce_sens depslatec)
target_link_libraries (pce_sens deplapack)
target_link_libraries (pce_sens depblas  )
target_link_libraries (pce_sens depfigtree )
target_link_libraries (pce_sens depann   )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pce_sens gfortran stdc++)
  else()
    target_link_libraries (pce_sens ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pce_sens ifcore)
  else()
    target_link_libraries (pce_sens ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_sens gfortran stdc++)
  else()
    target_link_libraries (pce_sens ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)

include_directories (${INSTALL_LIB_DIR})
include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
include_directories (/Users/kjohnst/UQTk-SNL/dep/sundials/include/ )
MESSAGE("{PROJECT_BINARY_DIR}/../../../dep/sundials/src = ${PROJECT_BINARY_DIR}/../../../dep/sundials/src")
#include_directories (${CMAKE_SUNDIALS_DIR}/include)
MESSAGE("{CMAKE_SUNDIALS_DIR}/include = ${CMAKE_SUNDIALS_DIR}/include")
include_directories(/Users/kjohnst/UQTk-build/dep/sundials/include)

INSTALL(TARGETS pce_sens DESTINATION bin)
