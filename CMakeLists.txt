cmake_minimum_required (VERSION 3.0)
project (UQTk)

if(CMAKE_VERSION VERSION_GREATER 3.14)
cmake_policy(SET CMP0086 OLD)
endif()

if(CMAKE_VERSION VERSION_GREATER 3.13)
cmake_policy(SET CMP0078 OLD)
endif()


# set( CMAKE_VERBOSE_MAKEFILE on ) # see all output


include( CTest )

# # set cpack for packaging
# SET(CPACK_GENERATOR "STGZ;TGZ;TZ")
# SET(CPACK_PACKAGE_NAME "UQTk")
# SET(CPACK_PACKAGE_VERSION "3.0")
# SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Bert Debusschere") #required
# INCLUDE(CPack)


IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX
    "${PROJECT_BINARY_DIR}" CACHE PATH "FOO install prefix" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(INSTALL_LIB_DIR     lib    )
set(INSTALL_BIN_DIR     bin    )
set(INSTALL_INCLUDE_DIR include)

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Add the lib directory for this package to the install RPATH
# (see https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling for more info)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Set it up so the installed binaries have the full RPATH so they can find all their libraries
# (see https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling for more info)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

option(PyUQTk "PyUQTk" OFF)
#option(DFI "DFI" OFF)

# CXX flags
#set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_CXX_FLAGS "-O2 -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HDF5 -I/opt/local/hdf5/include -L/opt/local/hdf5/lib -lhdf5_hl -lhdf5_cpp -lhdf5 -lhdf5_fortran -lhdf5_hl_cpp")

option(UQTk_ENABLE_SUNDIALS_SUBMODULE_UPDATE "Flag to enable sundials submodule update" OFF)

SET(BUILD_SUNDIALS OFF CACHE BOOL "BUILD_SUNDIALS is OFF by default")
if(NOT DEFINED CMAKE_SUNDIALS_DIR)
    UNSET(BUILD_SUNDIALS CACHE)
    SET(BUILD_SUNDIALS ON CACHE BOOL "To build sundials (vs specifying sundials directly)")
endif()

add_definitions(-D__wsu)
add_definitions(-fPIC)
add_definitions(-w)

add_subdirectory (dep     )
add_subdirectory (cpp     )
add_subdirectory (examples)
add_subdirectory (PyUQTk  )
add_subdirectory (doc  )

find_package(BLAS)
find_package(LAPACK)

if(NOT LAPACK_FOUND)
  MESSAGE("-- Not able to find system installed lapack. Using UQTk version of Lapack. If system version of Lapack is desired, please specify LAPACK_LIBRARIES in config file. ")
endif()

if(NOT BLAS_FOUND)
  MESSAGE("-- Not able to find system installed Blas. Using UQTk version of Blas. If system version of Blas is desired, please specify BLAS_LIBRARIES in config file. ")
endif()

# Copy over README file too
INSTALL(FILES README.md
        PERMISSIONS OWNER_WRITE OWNER_READ
        DESTINATION .)

# Copy over LICENSE file too
INSTALL(FILES LICENSE
        PERMISSIONS OWNER_WRITE OWNER_READ
        DESTINATION .)
