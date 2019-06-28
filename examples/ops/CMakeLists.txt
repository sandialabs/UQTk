project (UQTk)

SET(copy_FILES
  plot_logs.py
  plot_pdf.py
  LogComp.cpp
  Ops.cpp
  Utils.cpp
  Utils.h
  )


# Ops.x
add_executable (Ops.x Ops.cpp Utils.cpp)
add_executable (LogComp.x LogComp.cpp Utils.cpp)

target_link_libraries (Ops.x uqtkpce  )
target_link_libraries (Ops.x uqtkquad )
target_link_libraries (Ops.x uqtktools)
target_link_libraries (Ops.x uqtkarray )

target_link_libraries (Ops.x depdsfmt )
target_link_libraries (Ops.x depcvode )
target_link_libraries (Ops.x depnvec  )
target_link_libraries (Ops.x depslatec)
target_link_libraries (Ops.x deplapack)
target_link_libraries (Ops.x depblas  )
target_link_libraries (Ops.x depfigtree )
target_link_libraries (Ops.x depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (Ops.x gfortran stdc++)
  else()
    target_link_libraries (Ops.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (Ops.x ifcore)
  else()
    target_link_libraries (Ops.x ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (Ops.x gfortran stdc++)
  else()
    target_link_libraries (Ops.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

# LogComp.x
target_link_libraries (LogComp.x uqtkpce  )
target_link_libraries (LogComp.x uqtkquad )
target_link_libraries (LogComp.x uqtktools)
target_link_libraries (LogComp.x uqtkarray)

target_link_libraries (LogComp.x depdsfmt )
target_link_libraries (LogComp.x depcvode )
target_link_libraries (LogComp.x depnvec  )
target_link_libraries (LogComp.x depslatec)
target_link_libraries (LogComp.x deplapack)
target_link_libraries (LogComp.x depblas  )
target_link_libraries (LogComp.x depfigtree )
target_link_libraries (LogComp.x depann   )

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (LogComp.x gfortran stdc++)
  else()
    target_link_libraries (LogComp.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (LogComp.x ifcore)
  else()
    target_link_libraries (LogComp.x ${IntelLibPath}/libifcore.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (LogComp.x gfortran stdc++)
  else()
    target_link_libraries (LogComp.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../cpp/lib/include)
include_directories(../../cpp/lib/array)
include_directories(../../cpp/lib/pce     )
include_directories(../../cpp/lib/quad   )
include_directories(../../cpp/lib/tools  )

include_directories(../../dep/dsfmt)
include_directories(../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../dep/cvode-2.7.0/include")

INSTALL(TARGETS Ops.x LogComp.x 
              PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
              DESTINATION examples/ops)

INSTALL(FILES ${copy_FILES} 
              PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
              DESTINATION examples/ops)
