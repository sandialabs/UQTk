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

# # set cpack for packagin
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
set(INSTALL_INCLUDE_DIR include)

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

option(PyUQTk "PyUQTk" OFF)
#option(DFI "DFI" OFF)

# CXX flags
#set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_CXX_FLAGS "-O2 -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HDF5 -I/opt/local/hdf5/include -L/opt/local/hdf5/lib -lhdf5_hl -lhdf5_cpp -lhdf5 -lhdf5_fortran -lhdf5_hl_cpp")

add_definitions(-D__wsu)
add_definitions(-fPIC)
add_definitions(-w)

add_subdirectory (dep     )
add_subdirectory (cpp     )
add_subdirectory (examples)
add_subdirectory (PyUQTk  )
add_subdirectory (doc  )

# Copy over README file too
INSTALL(FILES README
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION .)

# Copy over LICENSE file too
INSTALL(FILES LICENSE
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION .)
