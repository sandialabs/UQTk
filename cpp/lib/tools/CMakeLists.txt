project(UQTk)

enable_language(Fortran)

SET(tools_HEADERS
  combin.h
  gq.h
  minmax.h
  multiindex.h
  pcmaps.h
  probability.h
  rosenblatt.h
  func.h
  tools.h
  )

FILE(GLOB slsrc "*.cpp")
FILE(GLOB fsrc  toolsf.f)

add_library(uqtktools ${fsrc} ${slsrc})

include_directories (../include)
include_directories (../array  )
include_directories (../pce)
include_directories (../quad)

include_directories (../../../dep/dsfmt )
include_directories (../../../dep/lapack)
include_directories (../../../dep/figtree)
include_directories (../../../dep/sundials/include )

include_directories ("${PROJECT_BINARY_DIR}/../../../dep/sundials/include")
include_directories (/Users/kjohnst/UQTk-SNL/dep/sundials/include/ )
MESSAGE("{PROJECT_BINARY_DIR}/../../../dep/sundials/src = ${PROJECT_BINARY_DIR}/../../../dep/sundials/src")
#include_directories (${CMAKE_SUNDIALS_DIR}/include)
MESSAGE("{CMAKE_SUNDIALS_DIR}/include = ${CMAKE_SUNDIALS_DIR}/include")

# Install the library
INSTALL(TARGETS uqtktools DESTINATION lib)

# Install the header files
INSTALL(FILES ${tools_HEADERS} DESTINATION include/uqtk)
