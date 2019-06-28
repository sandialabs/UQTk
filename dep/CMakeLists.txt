project (UQTk)

enable_language(Fortran)
enable_language(CXX)

add_subdirectory (slatec)
add_subdirectory (blas)
add_subdirectory (lapack)
add_subdirectory (dsfmt) 
add_subdirectory (lbfgs)
add_subdirectory (cvode-2.7.0)
add_subdirectory (ann)
add_subdirectory (figtree)

FILE(GLOB slatecsrc "slatec/*.f")
FILE(GLOB blassrc "blas/*.f")
FILE(GLOB lapacksrc "lapack/*.f")
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

add_library(depuqtk ${slatecsrc} ${blassrc} ${lapacksrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})

INSTALL(TARGETS depuqtk DESTINATION lib)
