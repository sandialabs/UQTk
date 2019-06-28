# ---------------------------------------------------------------
# $Revision: 1.3 $
# $Date: 2009/02/17 02:58:48 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# CMakeLists.txt file for the serial NVECTOR library

INSTALL(CODE "MESSAGE(\"\nInstall NVECTOR_SERIAL\n\")")

# Add variable nvecserial_SOURCES with the sources for the NVECSERIAL lib
SET(nvecserial_SOURCES nvector_serial.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECSERIAL library
SET(shared_SOURCES sundials_math.c)
ADD_PREFIX(${sundials_SOURCE_DIR}/sundials/ shared_SOURCES)

# Add variable nvecserial_HEADERS with the exported NVECSERIAL header files
SET(nvecserial_HEADERS nvector_serial.h)
ADD_PREFIX(${sundials_SOURCE_DIR}/include/nvector/ nvecserial_HEADERS)

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
ADD_LIBRARY(depnvec STATIC ${nvecserial_SOURCES} ${shared_SOURCES})
INSTALL(TARGETS depnvec DESTINATION lib)

# Install the NVECSERIAL header files
INSTALL(FILES ${nvecserial_HEADERS} DESTINATION include/dep/nvector)

# If FCMIX is enabled, build and install the FNVECSERIAL library
IF(FCMIX_ENABLE AND F77_FOUND)
  SET(fnvecserial_SOURCES fnvector_serial.c)
  IF(BUILD_STATIC_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_static STATIC ${fnvecserial_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecserial_static
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    INSTALL(TARGETS sundials_fnvecserial_static DESTINATION lib)
  ENDIF(BUILD_STATIC_LIBS)
  IF(BUILD_SHARED_LIBS)
    ADD_LIBRARY(sundials_fnvecserial_shared ${fnvecserial_SOURCES})
    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared
      PROPERTIES OUTPUT_NAME sundials_fnvecserial CLEAN_DIRECT_OUTPUT 1)
    SET_TARGET_PROPERTIES(sundials_fnvecserial_shared 
      PROPERTIES VERSION ${nveclib_VERSION} SOVERSION ${nveclib_SOVERSION})
    INSTALL(TARGETS sundials_fnvecserial_shared DESTINATION lib)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(FCMIX_ENABLE AND F77_FOUND)

#
MESSAGE(STATUS "Added NVECTOR_SERIAL module")
