project(UQTk)

SET(pce_HEADERS
  PCSet.h
  PCBasis.h
  )

add_library(uqtkpce PCBasis.cpp PCSet.cpp)

include_directories (../include)
include_directories (../array  )
include_directories (../tools  )
include_directories (../quad   )

# TPL
include_directories (../../../dep/slatec)
include_directories (../../../dep/dsfmt )
include_directories (../../../dep/sundials/include )
include_directories (${INSTALL_LIB_DIR})
include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
include_directories (/Users/kjohnst/UQTk-SNL/dep/sundials/include/ )
MESSAGE("{PROJECT_BINARY_DIR}/../../../dep/sundials/src = ${PROJECT_BINARY_DIR}/../../../dep/sundials/src")
#include_directories (${CMAKE_SUNDIALS_DIR}/include)
MESSAGE("{CMAKE_SUNDIALS_DIR}/include = ${CMAKE_SUNDIALS_DIR}/include")
include_directories(/Users/kjohnst/UQTk-build/dep/sundials/include)


# Install the library
INSTALL(TARGETS uqtkpce DESTINATION lib)

# Install the header files
INSTALL(FILES ${pce_HEADERS} DESTINATION include/uqtk)
