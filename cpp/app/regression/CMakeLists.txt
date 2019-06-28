project (UQTk)

add_executable (regression regression.cpp)

target_link_libraries (regression uqtklreg )
target_link_libraries (regression uqtkpce )
target_link_libraries (regression uqtkbcs  )
target_link_libraries (regression uqtkquad )
target_link_libraries (regression uqtktools)
target_link_libraries (regression uqtkarray)

target_link_libraries (regression depdsfmt )
target_link_libraries (regression deplbfgs )
target_link_libraries (regression depcvode )
target_link_libraries (regression depnvec  )
target_link_libraries (regression depslatec)
target_link_libraries (regression deplapack)
target_link_libraries (regression depblas  )
target_link_libraries (regression depfigtree )
target_link_libraries (regression depann   )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (regression gfortran stdc++)
  else()
    target_link_libraries (regression ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (regression ifcore)
    target_link_libraries (regression ifport)
  else()
    target_link_libraries (regression ${IntelLibPath}/libifcore.a)
    target_link_libraries (regression ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (regression gfortran stdc++)
  else()
    target_link_libraries (regression ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/bcs    )
include_directories(../../lib/lreg   )


include_directories(../../../dep/dsfmt)
include_directories(../../../dep/lbfgs)
include_directories(../../../dep/figtree)
include_directories(../../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS regression DESTINATION bin)

