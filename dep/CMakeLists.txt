project (UQTk)

enable_language(Fortran)
enable_language(CXX)

add_subdirectory (slatec)
add_subdirectory (dsfmt)
add_subdirectory (lbfgs)
# Sundials
if(NOT DEFINED CMAKE_SUNDIALS_DIR)

    FIND_PACKAGE(GIT)

    #set path to sundials
    SET(CMAKE_SUNDIALS_DIR ${CMAKE_INSTALL_PREFIX})

    SET(BUILD_SUNDIALS ON CACHE BOOL "To build sundials (vs specifying sundials directly)")

    # Process submodule Sundials (includes, cxx flags and lib)
    MESSAGE("-- UQTk uses submodule of sundials located at dep/sundials")
    # Automatically update sundials GIT
    # For development purpose inside of sundials, we can disable this update
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

    #turn off shared libraries
    SET(BUILD_SHARED_LIBS OFF CACHE BOOL "Build sundials shared library")
    #turn off example generation of all kinds
    SET(EXAMPLES_ENABLE_C OFF CACHE BOOL "Sundials examples in c")
    SET(EXAMPLES_ENABLE_F77 OFF CACHE BOOL "Sundials examples in fortran 77")
    SET(EXAMPLES_ENABLE_CXX OFF CACHE BOOL "Sundials examples in c++")
    SET(EXAMPLES_ENABLE_F90 OFF CACHE BOOL "Sundials examples in fortran 90")
    SET(EXAMPLES_ENABLE_CUDA OFF CACHE BOOL "Sundials examples in cuda")
    #build sundials
    ADD_SUBDIRECTORY(sundials)
    INCLUDE_DIRECTORIES(sundials/include)
    include_directories("${PROJECT_BINARY_DIR}/sundials/include")
endif()
add_subdirectory (ann)
add_subdirectory (figtree)

FILE(GLOB slatecsrc "slatec/*.f")
FILE(GLOB dsfmtsrc "dsfmt/*.c")
FILE(GLOB lbfgsCsrc "lbfgs/lbfgsDR.c")
FILE(GLOB lbfgsFsrc "lbfgs/lbfgs_routines.f")
FILE(GLOB annsrc "ann/*.cpp")
FILE(GLOB figtreesrc "figtree/*.cpp")


include_directories (dsfmt)
include_directories (lbfgs)
include_directories (ann)
include_directories (figtree)
include_directories("../cpp/lib/include")

find_package(BLAS)
find_package(LAPACK)

#if we need lapack or blas
if(NOT LAPACK_FOUND)
  add_subdirectory(lapack)
  FILE(GLOB lapacksrc "lapack/*.f")
endif()

if(NOT BLAS_FOUND)
  add_subdirectory(blas)
  FILE(GLOB blassrc "blas/*.f")
endif()


if(NOT BLAS_FOUND AND NOT LAPACK_FOUND)
  add_library(depuqtk ${lapacksrc} ${blassrc} ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
elseif(NOT BLAS_FOUND)
  add_library(depuqtk ${blassrc} ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
elseif(NOT LAPACK_FOUND)
  add_library(depuqtk ${lapacksrc} ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
else()
  add_library(depuqtk ${slatecsrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})
endif()



INSTALL(TARGETS depuqtk DESTINATION lib)
