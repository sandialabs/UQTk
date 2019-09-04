project (UQTk)

enable_language(Fortran)
enable_language(CXX)

add_subdirectory (slatec)
add_subdirectory (blas)
add_subdirectory (lapack)
add_subdirectory (dsfmt)
add_subdirectory (lbfgs)
# Sundials
if(NOT DEFINED CMAKE_SUNDIALS_DIR)

    FIND_PACKAGE(GIT)

    # Process submodule Sundials (includes, cxx flags and lib)
    MESSAGE("-- UQTk uses submodule of sundials located at dep/sundials")
    # Automatically update sundials GIT
    # For development purpose inside of sundials, we can disable this update
    MESSAGE("PROJECT_SOURCE_DIR = ${PROJECT_SOURCE_DIR}")
    IF (UQTk_ENABLE_SUNDIALS_SUBMODULE_UPDATE)
      IF (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        MESSAGE(STATUS "Git submodule updating")
        EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
                        RESULT_VARIABLE GIT_SUBMODULE_RESULT)
        IF (NOT GIT_SUBMODULE_RESULT EQUAL "0")
          MESSAGE(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMODULE_RESULT}")
        ENDIF()
      ELSE()
        MESSAGE(FATAL_ERROR "Git submodule updating fails ")
      ENDIF()
    ENDIF()

    #build sundials
    ADD_SUBDIRECTORY(dep/sundials)
    INCLUDE_DIRECTORIES(dep/sundials/include)
    MESSAGE("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    MESSAGE("sundials_INCLUDE_DIRS_RET = ${sundials_INCLUDE_DIRS_RET}")
    include_directories("${PROJECT_BINARY_DIR}/dep/sundials/include")
endif()
add_subdirectory (ann)
add_subdirectory (figtree)

FILE(GLOB slatecsrc "slatec/*.f")
FILE(GLOB blassrc "blas/*.f")
FILE(GLOB lapacksrc "lapack/*.f")
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

add_library(depuqtk ${slatecsrc} ${blassrc} ${lapacksrc} ${dsfmtsrc} ${lbfgsCsrc} ${lbfgsFsrc} ${annsrc} ${figtreesrc})

INSTALL(TARGETS depuqtk DESTINATION lib)
