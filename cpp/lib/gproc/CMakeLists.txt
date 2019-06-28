project(UQTk)

SET(gproc_HEADERS
  gproc.h
  )

add_library(uqtkgproc gproc.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )
include_directories (../pce    )
include_directories (../bcs    )


include_directories (../../../dep/lapack)
include_directories (../../../dep/blas)
include_directories (../../../dep/lbfgs)
include_directories (../../../dep/dsfmt)
include_directories (../../../dep/figtree)
include_directories (../../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../../dep/cvode-2.7.0/include")

# Install the library
INSTALL(TARGETS uqtkgproc DESTINATION lib)

# Install the header files
INSTALL(FILES ${gproc_HEADERS} DESTINATION include/uqtk)
