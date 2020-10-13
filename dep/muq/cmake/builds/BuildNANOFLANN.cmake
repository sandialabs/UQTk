set(EIGEN_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/nanoflann/src")

include(ExternalProject)
if( NOT NANOFLANN_EXTERNAL_SOURCE )
    set(NANOFLANN_EXTERNAL_SOURCE https://github.com/jlblancoc/nanoflann/archive/master.zip)
endif()

set(NANOFLANN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/muq_external/)
ExternalProject_Add(
 NANOFLANN
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/nanoflann
        URL ${NANOFLANN_EXTERNAL_SOURCE}
        BUILD_COMMAND ""
      	CONFIGURE_COMMAND ""
        INSTALL_COMMAND mkdir -p ${NANOFLANN_INSTALL_DIR}include && cp ${CMAKE_CURRENT_BINARY_DIR}/external/nanoflann/src/NANOFLANN/include/nanoflann.hpp ${NANOFLANN_INSTALL_DIR}include/
)

set_property( TARGET NANOFLANN PROPERTY FOLDER "Externals")

set(NANOFLANN_INCLUDE_DIRS ${NANOFLANN_INSTALL_DIR}/include/nanoflann)
message(STATUS "Adding ${NANOFLANN_INSTALL_DIR}/include/nanoflann for a nanoflann include directory.")
