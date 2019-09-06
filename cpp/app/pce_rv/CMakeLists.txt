project (UQTk)

link_directories(${CMAKE_SUNDIALS_DIR}/lib)
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/cvode")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/nvector/serial")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/sunlinsol/dense")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/sunmatrix/dense")

add_executable (pce_rv pce_rv.cpp)

target_link_libraries (pce_rv uqtkpce  )
target_link_libraries (pce_rv uqtkquad )
target_link_libraries (pce_rv uqtktools)
target_link_libraries (pce_rv uqtkarray)

target_link_libraries (pce_rv depdsfmt )
target_link_libraries (pce_rv sundials_cvode)
target_link_libraries (pce_rv sundials_nvecserial)
target_link_libraries (pce_rv sundials_sunlinsoldense)
target_link_libraries (pce_rv sundials_sunmatrixdense)
target_link_libraries (pce_rv depslatec)
target_link_libraries (pce_rv deplapack)
target_link_libraries (pce_rv depblas  )
target_link_libraries (pce_rv depfigtree )
target_link_libraries (pce_rv depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pce_rv gfortran stdc++)
  else()
    target_link_libraries (pce_rv ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pce_rv ifcore)
  else()
    target_link_libraries (pce_rv ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_rv gfortran stdc++)
  else()
    target_link_libraries (pce_rv ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS pce_rv DESTINATION bin)
