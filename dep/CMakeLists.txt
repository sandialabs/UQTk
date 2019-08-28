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

add_library(depuqtk ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})

INSTALL(TARGETS depuqtk DESTINATION lib)
