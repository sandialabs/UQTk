project(UQTk)

SET(lreg_HEADERS
  lreg.h
  )

add_library(uqtklreg lreg.cpp)

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
include_directories (../../../dep/sundials/include )
include_directories (${INSTALL_LIB_DIR})
include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
include_directories (/Users/kjohnst/UQTk-SNL/dep/sundials/include/ )
MESSAGE("{PROJECT_BINARY_DIR}/../../../dep/sundials/src = ${PROJECT_BINARY_DIR}/../../../dep/sundials/src")
#include_directories (${CMAKE_SUNDIALS_DIR}/include)
MESSAGE("{CMAKE_SUNDIALS_DIR}/include = ${CMAKE_SUNDIALS_DIR}/include")
include_directories(/Users/kjohnst/UQTk-build/dep/sundials/include)

# Install the library
INSTALL(TARGETS uqtklreg DESTINATION lib)

# Install the header files
INSTALL(FILES ${lreg_HEADERS} DESTINATION include/uqtk)
