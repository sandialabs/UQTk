include(ExternalProject)

set(GLOG_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/muq_external/)
if(DEFINED GLOG_EXTERNAL_SOURCE)
    ExternalProject_Add(
  	  GLOG
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/glog
                    URL ${GLOG_EXTERNAL_SOURCE}
                    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/external/glog/src/GLOG/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${CMAKE_CXX_FLAGS} --prefix=${GLOG_INSTALL_DIR}
                    BUILD_COMMAND make
                    BUILD_IN_SOURCE 1
                    INSTALL_COMMAND "make install"
    )
else()
    ExternalProject_Add(
  	  GLOG
                    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/external/glog
                    SVN_REPOSITORY http://google-glog.googlecode.com/svn/trunk/
                    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/external/glog/src/GLOG/configure CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CXXFLAGS=${CMAKE_CXX_FLAGS} --prefix=${GLOG_INSTALL_DIR}
                    BUILD_COMMAND make
                    BUILD_IN_SOURCE 1
                    INSTALL_COMMAND "make install"
    )
endif()


  set_property( TARGET GLOG PROPERTY FOLDER "Externals")
			   
  set( GLOG_INCLUDE_DIRS ${GLOG_INSTALL_DIR}/include )
  set( GLOG_LIBRARIES ${GLOG_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}glog${CMAKE_SHARED_LIBRARY_SUFFIX})
  set( GLOG_LIBRARIES_STATIC ${GLOG_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}glog${CMAKE_STATIC_LIBRARY_SUFFIX})

