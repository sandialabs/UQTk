project (UQTk)

add_executable (LowRankRegrTest main.cpp)

target_link_libraries (LowRankRegrTest uqtk )

target_link_libraries (LowRankRegrTest depdsfmt )
target_link_libraries (LowRankRegrTest m lapack ${LAPACK_LIBRARIES})
target_link_libraries (LowRankRegrTest m blas ${BLAS_LIBRARIES})
target_link_libraries (LowRankRegrTest depfigtree  )
target_link_libraries (LowRankRegrTest depann  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (LowRankRegrTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (LowRankRegrTest ifcore ifport)
  else()
    target_link_libraries (LowRankRegrTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (LowRankRegrTest ${IntelLibPath}/libifport.a)
  endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (LowRankRegrTest gfortran stdc++)
  else()
    target_link_libraries (LowRankRegrTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()


include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/tools  )
include_directories(../../lib/lowrank)

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/figtree)

INSTALL(TARGETS LowRankRegrTest DESTINATION bin/tests/)

add_test(LowRankRegrTest LowRankRegrTest)
