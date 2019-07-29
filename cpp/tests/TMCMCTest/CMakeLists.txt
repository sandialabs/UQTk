project (UQTk)

SET(copy_SCRIPT_FILES
  run.sh
  tmcmc_getLL.sh
  tmcmc_getLP.sh
  move.sh
  )

SET(copy_DATA_FILES
  rngs.dat
  )

SET(copy_EXE_FILES
    modelLLik.x
    modelLPr.x
   )

add_executable (modelLLik.x model.cpp)
add_executable (modelLPr.x lprior.cpp)
target_link_libraries (modelLPr.x depdsfmt)


INSTALL(FILES ${copy_FILES} DESTINATION cpp/tests/TMCMCTest/)

add_executable (TMCMCTest main.cpp)

target_link_libraries (TMCMCTest uqtk  )

target_link_libraries (TMCMCTest depdsfmt )
target_link_libraries (TMCMCTest depcvode )
target_link_libraries (TMCMCTest depnvec  )
target_link_libraries (TMCMCTest depslatec)
target_link_libraries (TMCMCTest deplapack)
target_link_libraries (TMCMCTest depblas  )
target_link_libraries (TMCMCTest deplbfgs  )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  target_link_libraries (TMCMCTest gfortran expat stdc++)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (TMCMCTest ifcore ifport)
  else()
    target_link_libraries (TMCMCTest ${IntelLibPath}/libifcore.a)
    target_link_libraries (TMCMCTest ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (TMCMCTest gfortran stdc++)
  else()
    target_link_libraries (TMCMCTest ${ClangLibPath}/libgfortran.dylib  ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

link_directories(${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/src/cvode)

include_directories(../../lib/pce    )
include_directories(../../lib/array  )
include_directories(../../lib/include)
include_directories(../../lib/quad   )
include_directories(../../lib/tools  )
include_directories(../../lib/tmcmc2  )

include_directories(../../../dep/dsfmt)
include_directories(../../../dep/lapack)
include_directories(../../../dep/blas)
include_directories(../../../dep/slatec)
include_directories(../../../dep/lbfgs)
include_directories(../../../dep/cvode-2.7.0/include)

include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

INSTALL(TARGETS TMCMCTest DESTINATION cpp/tests/TMCMCTest)

add_test(TMCMCTest TMCMCTest)

INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION cpp/tests/TMCMCTest
)

INSTALL(FILES ${copy_SCRIPT_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION cpp/tests/TMCMCTest
)

INSTALL(FILES ${copy_DATA_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION cpp/tests/TMCMCTest
)
