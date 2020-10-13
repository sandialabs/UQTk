include(ExternalProject)

if(NOT DEFINED BOOST_EXTERNAL_SOURCE)
  #NOTE: when updated to 1.63.0, you can remove the serialization patch below.
  set(BOOST_EXTERNAL_SOURCE http://downloads.sourceforge.net/project/boost/boost/1.63.0/boost_1_63_0.tar.gz)
endif()



set(BOOST_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/boost/src/BOOST")

# decide what toolset boost should use, start with compiler types, then work through operating systems
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

  # using Intel C++
  if(WIN32)
    set(BOOST_TOOLSET_NAME intel-win32)
  else()
    set(BOOST_TOOLSET_NAME intel-linux)
  endif()

  set(BOOST_CXX_FLAGS "-std=c++11")
  set(BOOST_LINK_FLAGS "")

# is this an OSX machine?
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  set(BOOST_TOOLSET_NAME darwin)

  if(MUQ_USE_LIBC11)
  	set(BOOST_CXX_FLAGS "-std=c++11 -stdlib=libc++")
  	set(BOOST_LINK_FLAGS "-stdlib=libc++")
  else(MUQ_USE_LIBC11)
	set(BOOST_CXX_FLAGS "-std=c++11")
	set(BOOST_LINK_FLAGS "")
  endif(MUQ_USE_LIBC11)

# is the compiler clang?
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

    set(BOOST_TOOLSET_NAME clang)

    if(MUQ_USE_LIBC11)
    	set(BOOST_CXX_FLAGS "-std=c++11 -stdlib=libc++")
    	set(BOOST_LINK_FLAGS "-stdlib=libc++")
    else(MUQ_USE_LIBC11)
        set(BOOST_CXX_FLAGS "-std=c++11")
        set(BOOST_LINK_FLAGS "")
    endif(MUQ_USE_LIBC11)

# is the compiler g++?
elseif(CMAKE_COMPILER_IS_GNUCXX)

    set(BOOST_TOOLSET_NAME gcc)

    set(BOOST_CXX_FLAGS "-std=c++11")
    set(BOOST_LINK_FLAGS "")

# is this an windows machine?
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

  if(MINGW)
    #using MINGW
    set(BOOST_TOOLSET_NAME mingw)
  elseif(MSYS)
    # using Visual Studio C++
    set(BOOST_TOOLSET_NAME msvc)
  else()
    message( FATAL_ERROR "Unable to find a BOOST toolset that matches your compiler and system.  Either use a different compiler, or try installing boost manually." )
  endif()

  set(BOOST_CXX_FLAGS "-std=c++11")
  set(BOOST_LINK_FLAGS "")

else()
        message( FATAL_ERROR "Unable to find a BOOST toolset that matches your compiler and system.  Either use a different compiler, or try installing boost manually." )
endif()


set(Boost_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/muq_external)

# create a configure file for boost_1_55_0
message(STATUS "Creating ${BOOST_BUILD_DIR}/tools/build/v2/user-config.jam")
string(REGEX MATCH "[0-9]+\\.[0-9]+" BOOST_TOOLSET_VERSION "${CMAKE_CXX_COMPILER_VERSION}")
if(MUQ_USE_OPENMPI)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/builds/user-config-mpi.jam.in ${CMAKE_CURRENT_BINARY_DIR}/user-config.jam)
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/builds/user-config.jam.in ${CMAKE_CURRENT_BINARY_DIR}/user-config.jam)
endif()

message(STATUS "BOOST_LINK_FLAGS = ${BOOST_LINK_FLAGS}")
message(STATUS "BOOST_CXX_FLAGS = ${BOOST_CXX_FLAGS}")


if(MUQ_USE_LIBC11)
    ExternalProject_Add(
        BOOST
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/boost
        URL ${BOOST_EXTERNAL_SOURCE}
        PATCH_COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/external/boost/shared_ptr_helper.hpp ${CMAKE_CURRENT_BINARY_DIR}/external/boost/src/BOOST/boost/serialization/shared_ptr_helper.hpp
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND mv ${CMAKE_CURRENT_BINARY_DIR}/user-config.jam ${BOOST_BUILD_DIR}/user-config.jam && ${BOOST_BUILD_DIR}/bootstrap.sh --prefix=${Boost_INSTALL_DIR} --without-icu
        BUILD_COMMAND ${BOOST_BUILD_DIR}/b2 cxxflags=${BOOST_CXX_FLAGS} linkflags=${BOOST_LINK_FLAGS} variant=release --user-config=${BOOST_BUILD_DIR}/user-config.jam toolset=${BOOST_TOOLSET_NAME}-muq  --with-filesystem --with-graph --with-system --disable-icu install
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
        )
else(MUQ_USE_LIBC11)
        ExternalProject_Add(
          BOOST
          PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/boost
          URL ${BOOST_EXTERNAL_SOURCE}
          PATCH_COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/external/boost/shared_ptr_helper.hpp ${CMAKE_CURRENT_BINARY_DIR}/external/boost/src/BOOST/boost/serialization/shared_ptr_helper.hpp
          UPDATE_COMMAND ""
          CONFIGURE_COMMAND mv ${CMAKE_CURRENT_BINARY_DIR}/user-config.jam ${BOOST_BUILD_DIR}/user-config.jam && ${BOOST_BUILD_DIR}/bootstrap.sh --prefix=${Boost_INSTALL_DIR} --without-icu
          BUILD_COMMAND ${BOOST_BUILD_DIR}/b2 cxxflags=${BOOST_CXX_FLAGS} variant=release --user-config=${BOOST_BUILD_DIR}/user-config.jam toolset=${BOOST_TOOLSET_NAME}-muq --with-filesystem --with-graph --disable-icu --with-system  install
          BUILD_IN_SOURCE 1
          INSTALL_COMMAND ""
          )
 endif(MUQ_USE_LIBC11)


set_property( TARGET BOOST PROPERTY FOLDER "Externals")

set( Boost_INCLUDE_DIRS ${Boost_INSTALL_DIR}/include )

set( Boost_LIBRARIES ${Boost_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}boost_filesystem${CMAKE_SHARED_LIBRARY_SUFFIX})
list(APPEND Boost_LIBRARIES ${Boost_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}boost_system${CMAKE_SHARED_LIBRARY_SUFFIX})
list(APPEND Boost_LIBRARIES ${Boost_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}boost_graph${CMAKE_SHARED_LIBRARY_SUFFIX})
list(APPEND Boost_LIBRARIES ${Boost_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}boost_regex${CMAKE_SHARED_LIBRARY_SUFFIX})

set( Boost_LIBRARIES_STATIC ${Boost_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}boost_filesystem${CMAKE_STATIC_LIBRARY_SUFFIX})
list(APPEND Boost_LIBRARIES_STATIC ${Boost_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}boost_system${CMAKE_STATIC_LIBRARY_SUFFIX})
list(APPEND Boost_LIBRARIES_STATIC ${Boost_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}boost_graph${CMAKE_STATIC_LIBRARY_SUFFIX})
list(APPEND Boost_LIBRARIES_STATIC ${Boost_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}boost_regex${CMAKE_STATIC_LIBRARY_SUFFIX})
