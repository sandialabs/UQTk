
find_package(PkgConfig)

if(NOT DEFINED MUQ_GTEST_DIR)
	pkg_check_modules(PC_GTEST QUIET libgtest)
	set(GTEST_DEFINITIONS ${PC_GTEST_CFLAGS_OTHER})

	find_path(GTEST_INCLUDE_DIR gtest/gtest.h
    		  HINTS ${PC_GTEST_INCLUDEDIR} ${PC_GTEST_INCLUDE_DIRS}
          	  PATH_SUFFIXES gtest )

	find_library(GTEST_LIBRARY NAMES gtest
           	     HINTS ${PC_GTEST_LIBDIR} ${PC_GTEST_LIBRARY_DIRS} )

	find_library(GTEST_LIBRARY_STATIC NAMES ${library_prefix}gtest.${static_library_suffix}
                HINTS ${PC_GTEST_LIBDIR} ${PC_GTEST_LIBRARY_DIRS} )
else()
	find_path(GTEST_INCLUDE_DIR gtest/gtest.h
	          HINTS ${MUQ_GTEST_DIR}/include
	          PATH_SUFFIXES gtest NO_DEFAULT_PATH)
			  
	find_library(GTEST_LIBRARY NAMES gtest
	             HINTS ${MUQ_GTEST_DIR}/lib ${MUQ_GTEST_DIR}/build NO_DEFAULT_PATH)

	find_library(GTEST_LIBRARY_STATIC NAMES ${library_prefix}gtest.${static_library_suffix}
	             HINTS ${MUQ_GTEST_DIR}/lib ${MUQ_GTEST_DIR}/build NO_DEFAULT_PATH)

endif()

set(GTEST_LIBRARIES_STATIC ${GTEST_LIBRARY_STATIC} )
			 
set(GTEST_LIBRARIES ${GTEST_LIBRARY} )
set(GTEST_INCLUDE_DIRS ${GTEST_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GTEST  DEFAULT_MSG
                                  GTEST_LIBRARY GTEST_INCLUDE_DIR)

mark_as_advanced(GTEST_INCLUDE_DIR GTEST_LIBRARY )

if( GTEST_LIBRARY OR GTEST_LIBRARY_STATIC )
    set(MUQ_HAS_GTEST 1)
else()
    set(MUQ_HAS_GTEST 0)
endif() 
