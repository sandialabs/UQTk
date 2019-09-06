project (UQTk)

link_directories(${CMAKE_SUNDIALS_DIR}/lib)
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/cvode")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/nvector/serial")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/sunlinsol/dense")
link_directories("${PROJECT_BINARY_DIR}/../../../dep/sundials/src/sunmatrix/dense")
add_executable (pdf_cl pdf_cl.cpp)

target_link_libraries (pdf_cl uqtktools)
target_link_libraries (pdf_cl uqtkarray)

target_link_libraries (pdf_cl depdsfmt )
target_link_libraries (pdf_cl sundials_cvode)
target_link_libraries (pdf_cl sundials_nvecserial)
target_link_libraries (pdf_cl sundials_sunlinsoldense)
target_link_libraries (pdf_cl sundials_sunmatrixdense)
target_link_libraries (pdf_cl depslatec)
target_link_libraries (pdf_cl deplapack)
target_link_libraries (pdf_cl depblas  )
target_link_libraries (pdf_cl depfigtree )
target_link_libraries (pdf_cl depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pdf_cl gfortran stdc++)
  else()
    target_link_libraries (pdf_cl ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pdf_cl ifcore)
  else()
    target_link_libraries (pdf_cl ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pdf_cl gfortran stdc++)
  else()
    target_link_libraries (pdf_cl ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

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

INSTALL(TARGETS pdf_cl DESTINATION bin)
