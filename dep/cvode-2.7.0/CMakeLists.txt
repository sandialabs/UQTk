# ---------------------------------------------------------------
# $Revision: 1.10 $
# $Date: 2010/12/15 22:28:16 $
# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# Copyright (c) 2007, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# ---------------------------------------------------------------
# Top level CMakeLists.txt for SUNDIALS (for cmake build system)


# -------------------------------------------------------------
# Initial commands
# -------------------------------------------------------------

# Require a fairly recent cmake version

CMAKE_MINIMUM_REQUIRED(VERSION 2.2)

# Project SUNDIALS (initially only C supported)

PROJECT(sundials C)

# Set some variables with info on the SUNDIALS project

SET(PACKAGE_BUGREPORT "radu@llnl.gov")
SET(PACKAGE_NAME "SUNDIALS")
SET(PACKAGE_STRING "SUNDIALS 2.4.0")
SET(PACKAGE_TARNAME "sundials")
SET(PACKAGE_VERSION "2.4.0")

# Prohibit in-source build

IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "In-source build prohibited.")
ENDIF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Hide/show some cache variables

MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(CLEAR
  CMAKE_C_COMPILER
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_DEBUG
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_RELWITHDEB)

# Specify the VERSION and SOVERSION for shared libraries

SET(cvodelib_VERSION "1.0.0")
SET(cvodelib_SOVERSION "1")

SET(nveclib_VERSION "0.0.2")
SET(nveclib_SOVERSION "0")

# Specify the location of additional CMAKE modules

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# -------------------------------------------------------------
# MACRO definitions
# -------------------------------------------------------------

# Macros to hide/show cached variables.
# These two macros can be used to "hide" or "show" in the
# list of cached variables various variables and/or options 
# that depend on other options.
# Note that once a variable is modified, it will preserve its
# value (hidding it merely makes it internal)

MACRO(HIDE_VARIABLE var)
  IF(DEFINED ${var})
    SET(${var} "${${var}}" CACHE INTERNAL "")
  ENDIF(DEFINED ${var})
ENDMACRO(HIDE_VARIABLE)

MACRO(SHOW_VARIABLE var type doc default)
  IF(DEFINED ${var})
    SET(${var} "${${var}}" CACHE "${type}" "${doc}" FORCE)
  ELSE(DEFINED ${var})
    SET(${var} "${default}" CACHE "${type}" "${doc}")
  ENDIF(DEFINED ${var})
ENDMACRO(SHOW_VARIABLE)

# Macros to append a common suffix or prefix to the elements of a list

MACRO(ADD_SUFFIX rootlist suffix)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${root}${suffix})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_SUFFIX)

MACRO(ADD_PREFIX prefix rootlist)
  SET(outlist )
  FOREACH(root ${${rootlist}})
    LIST(APPEND outlist ${prefix}${root})
  ENDFOREACH(root)
  SET(${rootlist} ${outlist})
ENDMACRO(ADD_PREFIX)

# Macro to print warning that some features will be disabled
# due to some failure.

MACRO(PRINT_WARNING message action)
  MESSAGE("WARNING: ${message}.\n   ${action}.")
ENDMACRO(PRINT_WARNING)

# Returns an unquoted string.  Note that CMake will readily turn such
# strings back into lists, due to the duality of lists and
# semicolon-separated strings.  So be careful how you use it.

MACRO(LIST2STRING alist astring)
  FOREACH(elem ${${alist}})
   SET(${astring} "${${astring}} ${elem}")
  ENDFOREACH(elem)
ENDMACRO(LIST2STRING)

# -------------------------------------------------------------
# Which modules to build?
# -------------------------------------------------------------
OPTION(BUILD_CVODE  "Build the CVODE library"  ON)


# -------------------------------------------------------------
# Other configuration options
# -------------------------------------------------------------

# Option that allows users to build static and/or shared libraries
# ----------------------------------------------------------------

OPTION(BUILD_STATIC_LIBS "Build static libraries" ON)
OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Prepare substitution variable SUNDIALS_EXPORT for sundials_config.h
# When building shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllexport)
# When linking to shared SUNDIALS libraries under Windows, use
#      #define SUNDIALS_EXPORT __declspec(dllimport)
# In all other cases (other platforms or static libraries
# under Windows), the SUNDIALS_EXPORT macro is empty

IF(BUILD_SHARED_LIBS AND WIN32)
  SET(SUNDIALS_EXPORT
    "#ifdef BUILD_SUNDIALS_LIBRARY
#define SUNDIALS_EXPORT __declspec(dllexport)
#else
#define SUNDIALS_EXPORT __declspec(dllimport)
#endif")
ELSE(BUILD_SHARED_LIBS AND WIN32)
  SET(SUNDIALS_EXPORT "#define SUNDIALS_EXPORT")
ENDIF(BUILD_SHARED_LIBS AND WIN32)

# Make sure we build at least one type of libraries
IF(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
  PRINT_WARNING("Both static and shared library generation were disabled"
    "Building static libraries was re-enabled")
  SET(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)
ENDIF(NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)

# Option to specify precision
# ---------------------------

SET(SUNDIALS_PRECISION "double" CACHE STRING "double, single or extended")

# prepare substitution variable PRECISION_LEVEL for sundials_config.h
STRING(TOUPPER ${SUNDIALS_PRECISION} SUNDIALS_PRECISION)
SET(PRECISION_LEVEL "#define SUNDIALS_${SUNDIALS_PRECISION}_PRECISION 1")

# Option to use the generic math libraries (UNIX only)
# ----------------------------------------------------

IF(UNIX)
  OPTION(USE_GENERIC_MATH "Use generic (std-c) math libraries" ON)
  IF(USE_GENERIC_MATH)
    # executables will be linked against -lm 
    SET(EXTRA_LINK_LIBS -lm)
    # prepare substitution variable GENERIC_MATH_LIB for sundials_config.h
    SET(GENERIC_MATH_LIB "#define SUNDIALS_USE_GENERIC_MATH")
  ENDIF(USE_GENERIC_MATH)
ENDIF(UNIX)

# -------------------------------------------------------------
# Enable Fortran support?
# -------------------------------------------------------------

# FCMIX support is an option only if at least one solver that 
# provides such an interface is built.

SHOW_VARIABLE(FCMIX_ENABLE BOOL "Enable Fortran-C support" OFF)

# -------------------------------------------------------------
# Enable BLAS/LAPACK support?
# -------------------------------------------------------------

OPTION(LAPACK_ENABLE "Enable Lapack support" OFF)

IF(NOT LAPACK_ENABLE)
  HIDE_VARIABLE(SUNDIALS_F77_FUNC_CASE)
  HIDE_VARIABLE(SUNDIALS_F77_FUNC_UNDERSCORES)
  HIDE_VARIABLE(LAPACK_LIBRARIES)
ENDIF(NOT LAPACK_ENABLE)

# -------------------------------------------------------------
# Enable MPI support?
# -------------------------------------------------------------

OPTION(MPI_ENABLE "Enable MPI support" OFF)

HIDE_VARIABLE(MPI_INCLUDE_PATH)
HIDE_VARIABLE(MPI_LIBRARIES)
HIDE_VARIABLE(MPI_EXTRA_LIBRARIES)
HIDE_VARIABLE(MPI_MPICC)
HIDE_VARIABLE(MPI_MPIF77)

# -------------------------------------------------------------
# Enable examples?
# -------------------------------------------------------------

OPTION(EXAMPLES_ENABLE "Build the SUNDIALS examples" OFF)
HIDE_VARIABLE(EXAMPLES_INSTALL)
HIDE_VARIABLE(EXAMPLES_INSTALL_PATH)
  
# -------------------------------------------------------------
# Add any other necessary compiler flags & definitions
# -------------------------------------------------------------

# -------------------------------------------------------------
# A Fortran compiler is needed if:
# (a) FCMIX is enabled
# (b) LAPACK is enabled (for the name-mangling scheme)
# -------------------------------------------------------------

IF(FCMIX_ENABLE OR LAPACK_ENABLE)
  INCLUDE(SundialsFortran)
  IF(NOT F77_FOUND AND FCMIX_ENABLE)
    PRINT_WARNING("Fortran compiler not functional"
      "FCMIX support will not be provided")
  ENDIF(NOT F77_FOUND AND FCMIX_ENABLE)
ENDIF(FCMIX_ENABLE OR LAPACK_ENABLE)

# -------------------------------------------------------------
# If we have a name-mangling scheme (either automatically
# inferred or provided by the user), set the SUNDIALS 
# compiler preprocessor macro definitions.
# -------------------------------------------------------------

SET(F77_MANGLE_MACRO1 "")
SET(F77_MANGLE_MACRO2 "")

IF(F77SCHEME_FOUND)
  # Symbols WITHOUT underscores
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub_")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub__")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) name ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "mysub__")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB_")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB__")
    SET(F77_MANGLE_MACRO1 "#define SUNDIALS_F77_FUNC(name,NAME) NAME ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MYSUB__")
  # Symbols with underscores 
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub_")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub__")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) name ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "my_sub__")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB_")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## _")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB_")
  IF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB__")
    SET(F77_MANGLE_MACRO2 "#define SUNDIALS_F77_FUNC_(name,NAME) NAME ## __")
  ENDIF(${CMAKE_Fortran_SCHEME_NO_UNDERSCORES} MATCHES "MY_SUB__")
ENDIF(F77SCHEME_FOUND)

# -------------------------------------------------------------
# Find (and test) the Lapack libraries
# -------------------------------------------------------------

# If LAPACK is needed, first try to find the appropriate
# libraries and linker flags needed to link against them.

# Macro to be inserted in sundials_config.h
SET(BLAS_LAPACK_MACRO "#define SUNDIALS_BLAS_LAPACK 0")

IF(LAPACK_ENABLE)

  INCLUDE(SundialsLapack)

  IF(LAPACK_FOUND)
    SET(BLAS_LAPACK_MACRO "#define SUNDIALS_BLAS_LAPACK 1")
  ELSE(LAPACK_FOUND)
    SHOW_VARIABLE(LAPACK_LIBRARIES STRING "Lapack libraries" "${LAPACK_LIBRARIES}")
    SHOW_VARIABLE(LAPACK_LINKER_FLAGS STRING "Lapack required linker flags" "${LAPACK_LINKER_FLAGS}")
  ENDIF(LAPACK_FOUND)

  IF(LAPACK_LIBRARIES AND NOT LAPACK_FOUND)
    PRINT_WARNING("LAPACK not functional"
      "Blas/Lapack support will not be provided")
  ENDIF(LAPACK_LIBRARIES AND NOT LAPACK_FOUND)

ELSE(LAPACK_ENABLE)

  HIDE_VARIABLE(LAPACK_LIBRARIES)
  HIDE_VARIABLE(LAPACK_LINKER_FLAGS)

ENDIF(LAPACK_ENABLE)

# -------------------------------------------------------------
# Configure the header file sundials_config.h
# -------------------------------------------------------------

# All required substitution variables should be available at this point.
# Generate the header file and place it in the binary dir.
CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_config.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  )

# Add the include directory in the source tree and the one in
# the binary tree (for the header file sundials_config.h)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)

# -------------------------------------------------------------
# Add selected modules to the build system
# -------------------------------------------------------------

# Shared components

ADD_SUBDIRECTORY(sundials)
ADD_SUBDIRECTORY(nvec_ser)
ADD_SUBDIRECTORY(cvode)

#----------------------------------
# Install configuration header file
#----------------------------------

# install configured header file
INSTALL(
  FILES ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  DESTINATION include/dep/sundials
  )

