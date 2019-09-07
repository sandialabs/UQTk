project (UQTk)

enable_language(Fortran)
enable_language(CXX)

add_subdirectory (slatec)
add_subdirectory (dsfmt)
add_subdirectory (lbfgs)
add_subdirectory (cvode-2.7.0)
add_subdirectory (ann)
add_subdirectory (figtree)

FILE(GLOB slatecsrc "slatec/*.f")
FILE(GLOB dsfmtsrc "dsfmt/*.c")
FILE(GLOB lbfgsCsrc "lbfgs/lbfgsDR.c")
FILE(GLOB lbfgsFsrc "lbfgs/lbfgs_routines.f")
FILE(GLOB annsrc "ann/*.cpp")
FILE(GLOB figtreesrc "figtree/*.cpp")


include_directories (dsfmt)
include_directories (lbfgs)
include_directories (ann)
include_directories (figtree)
include_directories("../cpp/lib/include")

find_package(BLAS)
find_package(LAPACK)


#if we need lapack or blas
if(NOT LAPACK_FOUND)
  add_subdirectory(lapack)
  FILE(GLOB lapacksrc "lapack/*.f")
endif()

if(NOT BLAS_FOUND)
  add_subdirectory(blas)
  FILE(GLOB blassrc "blas/*.f")
endif()


if(NOT BLAS_FOUND AND NOT LAPACK_FOUND)
  add_library(depuqtk ${lapacksrc} ${blassrc} ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
elseif(NOT BLAS_FOUND)
  add_library(depuqtk ${blassrc} ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
elseif(NOT LAPACK_FOUND)
  add_library(depuqtk ${lapacksrc} ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
else()
  add_library(depuqtk ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
endif()



INSTALL(TARGETS depuqtk DESTINATION lib)
