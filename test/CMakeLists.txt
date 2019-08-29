cmake_minimum_required (VERSION 3.0)
project (HELLO)


# set( CMAKE_VERBOSE_MAKEFILE on ) # see all output
include( CTest )

#Installation

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

option(PyUQTk "PyUQTk" OFF)
option(DFI "DFI" OFF)

# CXX flags
#set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_CXX_FLAGS "-O2 -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HDF5 -I/opt/local/hdf5/include -L/opt/local/hdf5/lib -lhdf5_hl -lhdf5_cpp -lhdf5 -lhdf5_fortran -lhdf5_hl_cpp")


#
# Options : use TCHEM prefix
#
OPTION(UQTk_ENABLE_CVODE "Flag to enable Cvode" ON)
OPTION(UQTk_ENABLE_CVODE_SUBMODULE_UPDATE "Flag to enable CVode submodule update" ON)
OPTION(CVODE_INSTALL_PATH "Path to Cvode installation")

#SET(HAVE_KOKKOS ${TCHEM_ENABLE_KOKKOS})

#
# Git
#
FIND_PACKAGE(Git)

#
#CVode
#

IF(CVODE_INSTALL_PATH)
  MESSAGE("-- Test uses CVode installed at ${CVODE_INSTALL_PATH}")
  include_directories(${CVODE_INSTALL_PATH})
ELSE()
  MESSAGE(" -- Please specify CVODE path. CVODE will not be linked")

  # Process submodule Kokkos (includes, cxx flags and lib)
  MESSAGE("-- TChem uses submodule of Kokkos located at external/sundials")
  # Automatically update Kokkos GIT
  # For development purpose inside of Kokkos, we can disable this update
  IF (UQTk_ENABLE_CVODE_SUBMODULE_UPDATE)
    IF (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/../.git")
      MESSAGE(STATUS "Git submodule updating")
      EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../
                      RESULT_VARIABLE GIT_SUBMODULE_RESULT)
      IF (NOT GIT_SUBMODULE_RESULT EQUAL "0")
        MESSAGE(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMODULE_RESULT}")
      ENDIF()
    ELSE()
      MESSAGE(FATAL_ERROR "Git submodule updating fails ")
    ENDIF()
  ENDIF()
  ADD_SUBDIRECTORY(external/sundials)
  MESSAGE("sundials_INCLUDE_DIRS_RET = ${sundials_INCLUDE_DIRS_RET}")
  INCLUDE_DIRECTORIES(external/sundials/include)
  MESSAGE("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
  include_directories("${PROJECT_BINARY_DIR}/external/sundials/include")

ENDIF()





add_executable(hello hello.cpp)
