project (UQTk)

SET(copy_SCRIPT_FILES
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  move.sh
  )

SET(copy_EXE_FILES
    model.x
   )

add_executable (model.x model.cpp)
add_executable (TMCMC2dTest main.cpp)

target_link_libraries (TMCMC2dTest uqtk  )

target_link_libraries (TMCMC2dTest depdsfmt )
target_link_libraries (TMCMC2dTest depcvode )
target_link_libraries (TMCMC2dTest depnvec  )
target_link_libraries (TMCMC2dTest depslatec)
target_link_libraries (TMCMC2dTest deplapack)
target_link_libraries (TMCMC2dTest depblas  )
target_link_libraries (TMCMC2dTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (TMCMC2dTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (TMCMC2dTest ifcore ifport)
  else()
    target_link_libraries (TMCMC2dTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (TMCMC2dTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (TMCMC2dTest gfortran stdc++)
  else()
    target_link_libraries (TMCMC2dTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
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

INSTALL(TARGETS TMCMC2dTest DESTINATION bin/tests/)

add_test(TMCMC2dTest TMCMC2dTest)


INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin/tests/
)

INSTALL(FILES ${copy_SCRIPT_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin/tests/
)

