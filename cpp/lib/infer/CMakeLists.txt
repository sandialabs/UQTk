project (UQTk)

SET(infer_HEADERS
  inference.h
  mrv.h
  post.h
  )

add_library(uqtkinfer inference.cpp mrv.cpp post.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )
include_directories (../pce    )
include_directories (../mcmc   )

include_directories (../../../dep/slatec)
include_directories (../../../dep/lapack)
include_directories (../../../dep/dsfmt )
include_directories (../../../dep/figtree )
include_directories (../../../dep/sundials/include )
include_directories (${INSTALL_LIB_DIR})
include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
include_directories (/Users/kjohnst/UQTk-SNL/dep/sundials/include/ )
MESSAGE("{PROJECT_BINARY_DIR}/../../../dep/sundials/src = ${PROJECT_BINARY_DIR}/../../../dep/sundials/src")
#include_directories (${CMAKE_SUNDIALS_DIR}/include)
MESSAGE("{CMAKE_SUNDIALS_DIR}/include = ${CMAKE_SUNDIALS_DIR}/include")


# Install the library
INSTALL(TARGETS uqtkinfer DESTINATION lib)

# Install the header files
INSTALL(FILES ${infer_HEADERS} DESTINATION include/uqtk)
