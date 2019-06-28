project (UQTk)

SET(copy_FILES
  pce_bcs.cpp
  pce_bcs.h
  )

# pce_bcs.x
add_executable (pce_bcs.x pce_bcs.cpp)

target_link_libraries (pce_bcs.x uqtkbcs  )
target_link_libraries (pce_bcs.x uqtkpce  )
target_link_libraries (pce_bcs.x uqtkmcmc )
target_link_libraries (pce_bcs.x uqtkquad )
target_link_libraries (pce_bcs.x uqtktools)
target_link_libraries (pce_bcs.x uqtkarray )
target_link_libraries (pce_bcs.x uqtkxmlutils )

target_link_libraries (pce_bcs.x depdsfmt )
target_link_libraries (pce_bcs.x deplbfgs )
target_link_libraries (pce_bcs.x depcvode )
target_link_libraries (pce_bcs.x depnvec  )
target_link_libraries (pce_bcs.x depslatec)
target_link_libraries (pce_bcs.x deplapack)
target_link_libraries (pce_bcs.x depblas  )
target_link_libraries (pce_bcs.x depfigtree )
target_link_libraries (pce_bcs.x depann   )

target_link_libraries (pce_bcs.x expat )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (pce_bcs.x gfortran stdc++)
  else()
    target_link_libraries (pce_bcs.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (pce_bcs.x ifcore ifport)
  else()
    target_link_libraries (pce_bcs.x ${IntelLibPath}/libifcore.a)
    target_link_libraries (pce_bcs.x ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (pce_bcs.x gfortran stdc++)
  else()
    target_link_libraries (pce_bcs.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/pce_bcs)
include_directories(../../cpp/lib/include    )
include_directories(../../cpp/lib/array      )
include_directories(../../cpp/lib/pce        )
include_directories(../../cpp/lib/bcs        )
include_directories(../../cpp/lib/quad       )
include_directories(../../cpp/lib/mcmc       )
include_directories(../../cpp/lib/tools      )
include_directories(../../cpp/lib/xmlutils   )

include_directories(../../dep/dsfmt)
include_directories(../../dep/figtree)
include_directories(../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../dep/cvode-2.7.0/include")

INSTALL(TARGETS pce_bcs.x 
              PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
              DESTINATION examples/pce_bcs)

INSTALL(FILES ${copy_FILES} DESTINATION examples/pce_bcs)
