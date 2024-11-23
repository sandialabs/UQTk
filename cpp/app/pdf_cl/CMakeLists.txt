project (UQTk)

add_executable (pdf_cl pdf_cl.cpp)

target_link_libraries (pdf_cl uqtktools)
target_link_libraries (pdf_cl uqtkarray)

target_link_libraries (pdf_cl depdsfmt )
target_link_libraries (pdf_cl sundials_core)
target_link_libraries (pdf_cl sundials_cvode)
target_link_libraries (pdf_cl sundials_nvecserial)
target_link_libraries (pdf_cl sundials_sunlinsoldense)
target_link_libraries (pdf_cl sundials_sunmatrixdense)
target_link_libraries (pdf_cl depslatec)
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
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pdf_cl gfortran stdc++)
  else()
    target_link_libraries (pdf_cl ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

target_link_libraries (pdf_cl m lapack ${LAPACK_LIBRARIES})
target_link_libraries (pdf_cl m blas ${BLAS_LIBRARIES})

include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)
include_directories (${CMAKE_SUNDIALS_DIR}/include)
if( BUILD_SUNDIALS)
  include_directories (../../../dep/sundials/include )
  include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
endif()

INSTALL(TARGETS pdf_cl DESTINATION bin)
