

include(CheckCXXCompilerFlag)

#use the MPI wrapper as necessary
option(MUQ_USE_MPI "Whether or not to compile MUQ with MPI support." OFF)
if(MUQ_USE_MPI)

  find_package(MPI REQUIRED)

  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_directories(${MPI_CXX_LIBRARIES})

  set(MUQ_HAS_MPI 1)

else(MUQ_USE_MPI)
  set(MUQ_HAS_MPI 0)
endif(MUQ_USE_MPI)

set(CMAKE_CXX_FLAGS_DEBUG  "-O0") #-O0 works better for memcheck
set(CMAKE_CXX_FLAGS_RELEASE  "-O3") #full optimization with debug symbols for profiling

set(CMAKE_CXX_FLAGS "-g")

# default to a release build
message(STATUS "User defined build type = " ${CMAKE_BUILD_TYPE})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()
message(STATUS "Final build type = " ${CMAKE_BUILD_TYPE})

set(MUQ_USE_LIBC11 OFF) # will turn on below if using clang and found

# check for c++11 support and add the required compiler flags
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

   execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)

   # check the gcc version to make sure it supports c++11
   if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        message(STATUS "C++11 found.")
   else ()
        message(FATAL_ERROR "A full implementation of C++11 is needed. When using g++ this means the gcc compiler must be 4.7 or newer.")
   endif()

   # check to make sure c++11 flag works
   CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_CXX11)
   if(NOT HAS_CXX11)
	   message(FATAL_ERROR "A check of the '-std=c++11' compiler flag flagged.  It seems that the compiler does not support c++11.")
   endif()

   # set compiler flags for g++
   set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -g -Wno-maybe-uninitialized -Wno-sign-compare -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-local-typedefs")


elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

    # check to make sure c++11 flag works
    CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_CXX11)
    if(NOT HAS_CXX11)
 	   message(FATAL_ERROR "A check of the '-std=c++11' compiler flag flagged.  It seems that the compiler does not support c++11.")
    endif()
	set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")

    CHECK_CXX_COMPILER_FLAG("-std=c++11 -stdlib=libc++" HAS_LIBCXX11)
    INCLUDE(LibcxxCheck)

    # set compiler flags for clang
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -g -Wno-unused-function -Wno-redeclared-class-member -Wno-deprecated-register -Wno-uninitialized -Wno-sign-compare -Wno-unknown-pragmas -Wunused-function -Wno-unused-variable -Wno-overloaded-virtual")

else()
    message(FATAL_ERROR "Your C++ compiler is not recognized or does not seem to support C++11.\nIf cmake did not find the correct compiler, try setting CMAKE_CXX_COMPILER to a suitable compiler.\n")

endif()

# Check to see if const& and by value need to be treated separately in AnyConstCast
INCLUDE(AnyCastCheck)

IF(MUQ_USE_OPENMP)
	CHECK_CXX_COMPILER_FLAG("-fopenmp" HAS_FOPENMP)
	CHECK_CXX_COMPILER_FLAG("-pthread" HAS_PTHREAD)

	if(HAS_FOPENMP AND HAS_PTHREAD)
  	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -pthread -ldl")
   else()
	   if(HAS_FOPENMP)
	     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -ldl")
	  else()
		  message(WARNING "The flag MUQ_USE_OPENMP is ON, but the compiler does not seem to support the -fopenmp flag.  OPENMP will not be used.")
	  endif()
  endif()
ENDIF(MUQ_USE_OPENMP)


# this is required for cmake version 3.0.0 and later
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
endif()

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
