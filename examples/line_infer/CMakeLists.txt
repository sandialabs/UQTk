project (UQTk)

SET(copy_FILES
  line_infer.cpp
  model.h
  posterior.h
  XMLreader.h
  line_infer.xml.templ
  line_infer.py
  )

# line_infer.x
add_executable (line_infer.x line_infer.cpp)

target_link_libraries (line_infer.x uqtkpce  )
target_link_libraries (line_infer.x uqtkmcmc )
target_link_libraries (line_infer.x uqtkquad )
target_link_libraries (line_infer.x uqtktools)
target_link_libraries (line_infer.x uqtkarray )
target_link_libraries (line_infer.x uqtkxmlutils )

target_link_libraries (line_infer.x depdsfmt )
target_link_libraries (line_infer.x deplbfgs )
target_link_libraries (line_infer.x depcvode )
target_link_libraries (line_infer.x depnvec  )
target_link_libraries (line_infer.x depslatec)
target_link_libraries (line_infer.x deplapack)
target_link_libraries (line_infer.x depblas  )
target_link_libraries (line_infer.x depfigtree  )
target_link_libraries (line_infer.x depann   )

target_link_libraries (line_infer.x expat )


# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (line_infer.x gfortran stdc++)
  else()
    target_link_libraries (line_infer.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (line_infer.x ifcore ifport)
  else()
    target_link_libraries (line_infer.x ${IntelLibPath}/libifcore.a)
    target_link_libraries (line_infer.x ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (line_infer.x gfortran stdc++)
  else()
    target_link_libraries (line_infer.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/line_infer)
include_directories(../../cpp/lib/include    )
include_directories(../../cpp/lib/array      )
include_directories(../../cpp/lib/pce        )
include_directories(../../cpp/lib/quad       )
include_directories(../../cpp/lib/mcmc       )
include_directories(../../cpp/lib/tools      )
include_directories(../../cpp/lib/xmlutils   )

include_directories(../../dep/dsfmt)
include_directories(../../dep/figtree)
include_directories(../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../dep/cvode-2.7.0/include")

INSTALL(TARGETS line_infer.x DESTINATION examples/line_infer)

INSTALL(FILES ${copy_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
        DESTINATION examples/line_infer)
