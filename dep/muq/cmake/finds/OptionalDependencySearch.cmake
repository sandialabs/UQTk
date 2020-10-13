########################################
##### LOOK FOR GTEST              ######
########################################
IF(MUQ_USE_GTEST)

  find_package(GTEST)
  set(MUQ_NEEDS_GTEST ON)

  IF(GTEST_FOUND)
    include(CheckGTEST)
  ENDIF()

  IF(GTEST_FOUND AND NOT GTEST_TEST_FAIL)
    set(MUQ_BUILD_TESTS ON)
    include_directories(${GTEST_INCLUDE_DIRS})
    LIST(APPEND test_link_libs ${GTEST_LIBRARIES})

  ELSE()

    message(WARNING "Could not find GTEST.  No tests can be compiled!")
    set(MUQ_BUILD_TESTS OFF)

  ENDIF(GTEST_FOUND AND NOT GTEST_TEST_FAIL)

ELSE(MUQ_USE_GTEST)

    message(WARNING “Tried to compile tests, but MUQ_USE_GTEST is OFF.  Turning off tests.”)
    set(MUQ_BUILD_TESTS OFF)
    set(MUQ_NEEDS_GTEST OFF)
ENDIF(MUQ_USE_GTEST)


########################################
##### LOOK FOR MKL                ######
########################################
if(MUQ_USE_MKL)

    # include the mkl library for linking
    if(MUQ_USE_OPENMP)
	  LIST(APPEND MUQ_LINK_LIBS ${MUQ_MKL_DIR}/lib/intel64/libmkl_intel_lp64${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_core${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_gnu_thread${CMAKE_SHARED_LIBRARY_SUFFIX})
    else()
	  LIST(APPEND MUQ_LINK_LIBS ${MUQ_MKL_DIR}/lib/intel64/libmkl_intel_lp64${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_core${CMAKE_SHARED_LIBRARY_SUFFIX} ${MUQ_MKL_DIR}/lib/intel64/libmkl_sequential${CMAKE_SHARED_LIBRARY_SUFFIX})
    endif()

    include_directories(${MUQ_MKL_DIR}/include)
    LIST(APPEND MUQ_EXTERNAL_INCLUDES ${MUQ_MKL_DIR}/include)
    add_definitions(-DEIGEN_USE_MKL_ALL)

endif()

########################################
##### LOOK FOR PYTHON             ######
########################################
list (FIND MUQ_REQUIRES PYTHON dindex)
if (${dindex} GREATER -1)
    set(MUQ_NEEDS_PYTHON ON)

    if(MUQ_USE_PYTHON)
        set(PYBIND11_CPP_STANDARD -std=c++11)

        FIND_PACKAGE(pybind11)

        if(NOT pybind11_FOUND)
            add_subdirectory(${CMAKE_SOURCE_DIR}/external/pybind11)
            include_directories(${CMAKE_SOURCE_DIR}/external/pybind11/include)
        endif()
    endif()
else()
    set(MUQ_NEEDS_PYTHON OFF)
    set(MUQ_USE_PYTHON OFF)
endif()

########################################
##### LOOK FOR DOLFIN/Fenics      ######
########################################
list (FIND MUQ_REQUIRES DOLFIN dindex)
message(${MUQ_REQUIRES})
if (${dindex} GREATER -1)
    set(MUQ_NEEDS_DOLFIN ON)

    if(MUQ_USE_DOLFIN AND NOT MUQ_USE_PYTHON)
        message(WARNING "Requested compilation with Fenics/Dolfin, but not Python.  Building the Fenics/Dolfin bindings requires building MUQ with Python support.")
        set(MUQ_USE_DOLFIN OFF)
    endif()

    if(MUQ_USE_DOLFIN)

        find_package(DOLFIN)

        if(DOLFIN_FOUND)
            if (EXISTS ${DOLFIN_USE_FILE})
                include(${DOLFIN_USE_FILE})

                # Default build type (can be overridden by user)
                if (NOT CMAKE_BUILD_TYPE)
                    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
                        "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
                endif()
            else()
                # Compiler definitions
                add_definitions(${DOLFIN_CXX_DEFINITIONS})

                # Include directories
                include_directories(${DOLFIN_INCLUDE_DIRS})
                include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
            endif()
        else()
            set(MUQ_USE_DOLFIN OFF)
        endif()

    endif()
else()
    set(MUQ_NEEDS_DOLFIN OFF)
    set(MUQ_USE_DOLFIN OFF)
endif()

########################################
##### REMOVE DUPLICATE INCLUDES   ######
########################################
if(MUQ_LINK_LIBS)
  list( REMOVE_DUPLICATES MUQ_LINK_LIBS)
endif()
