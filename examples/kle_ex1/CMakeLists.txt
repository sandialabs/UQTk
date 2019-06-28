project (UQTk)

SET(copy_SRC_FILES
    kl_1D.cpp
    kl_2D.cpp
    kl_2Du.cpp
    kl_utils.cpp
    kl_utils.h
   )

SET(copy_DAT_FILES
    data/cali_grid.dat
    data/cali_tria.dat
   )

SET(copy_SCRIPT_FILES
    run1D.sh
    run2D.sh
    mkplots.py
    pyutils.py
   )

SET(copy_EXE_FILES
    kl_1D.x
    kl_2D.x
    kl_2Du.x
   )

# kl_1D.x
add_executable (kl_1D.x kl_1D.cpp kl_utils.cpp)
target_link_libraries (kl_1D.x uqtktools uqtkarray uqtkkle)
target_link_libraries (kl_1D.x depdsfmt depslatec deplapack depblas depann depfigtree)

# kl_2D.x
add_executable (kl_2D.x kl_2D.cpp kl_utils.cpp)
target_link_libraries (kl_2D.x uqtktools uqtkarray uqtkkle)
target_link_libraries (kl_2D.x depdsfmt depslatec deplapack depblas depann depfigtree)

# kl_2Du.x
add_executable (kl_2Du.x kl_2Du.cpp kl_utils.cpp)
target_link_libraries (kl_2Du.x uqtktools uqtkarray uqtkkle)
target_link_libraries (kl_2Du.x depdsfmt depslatec deplapack depblas depann depfigtree)

# Link fortran libraries
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if ("${GnuLibPath}" STREQUAL "")
    target_link_libraries (kl_1D.x  gfortran stdc++)
    target_link_libraries (kl_2D.x  gfortran stdc++)
    target_link_libraries (kl_2Du.x gfortran stdc++)
  else()
    target_link_libraries (kl_1D.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
    target_link_libraries (kl_2D.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
    target_link_libraries (kl_2Du.x ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel
  if ("${IntelLibPath}" STREQUAL "")
    target_link_libraries (kl_1D.x    ifcore ifport)
    target_link_libraries (kl_2D.x    ifcore ifport)
    target_link_libraries (kl_2Du.x   ifcore ifport)
  else()
    target_link_libraries (kl_1D.x    ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
    target_link_libraries (kl_2D.x    ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
    target_link_libraries (kl_2Du.x   ${IntelLibPath}/libifcore.a ${IntelLibPath}/libifport.a)
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  if ("${ClangLibPath}" STREQUAL "")
    target_link_libraries (kl_1D.x  gfortran stdc++)
    target_link_libraries (kl_2D.x  gfortran stdc++)
    target_link_libraries (kl_2Du.x gfortran stdc++)
  else()
    target_link_libraries (kl_1D.x  ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
    target_link_libraries (kl_2D.x  ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
    target_link_libraries (kl_2Du.x ${ClangLibPath}/libgfortran.dylib ${ClangLibPath}/libquadmath.dylib ${ClangLibPath}/libstdc++.dylib)
  endif()
endif()

include_directories(../../examples/kle_ex1)
include_directories(../../cpp/lib/include )
include_directories(../../cpp/lib/array   )
include_directories(../../cpp/lib/tools   )
include_directories(../../cpp/lib/kle     )

include_directories(../../dep/dsfmt  )
include_directories(../../dep/lapack )
include_directories(../../dep/figtree)

INSTALL(TARGETS ${copy_EXE_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/kle_ex1
)
INSTALL(FILES ${copy_SRC_FILES}    
        PERMISSIONS OWNER_WRITE OWNER_READ
        DESTINATION examples/kle_ex1
)
INSTALL(FILES ${copy_SCRIPT_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/kle_ex1
)
INSTALL(FILES ${copy_DAT_FILES}    
        PERMISSIONS OWNER_WRITE OWNER_READ
        DESTINATION examples/kle_ex1/data
)

