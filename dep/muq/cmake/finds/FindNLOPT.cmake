
find_package(PkgConfig)

if(NOT DEFINED MUQ_NLOPT_DIR)
	pkg_check_modules(PC_NLOPT QUIET libnlopt)
	set(NLOPT_DEFINITIONS ${PC_NLOPT_CFLAGS_OTHER})

	find_path(NLOPT_INCLUDE_DIR nlopt.h
          HINTS ${PC_NLOPT_INCLUDEDIR} ${PC_NLOPT_INCLUDE_DIRS}
          PATH_SUFFIXES nlopt )

	find_library(NLOPT_LIBRARY NAMES nlopt nlopt_cxx
             HINTS ${PC_NLOPT_LIBDIR} ${PC_NLOPT_LIBRARY_DIRS} )

	find_library(NLOPT_LIBRARY_STATIC NAMES ${library_prefix}nlopt.${static_library_suffix} ${library_prefix}nlopt_cxx.${static_library_suffix}
             HINTS ${PC_NLOPT_LIBDIR} ${PC_NLOPT_LIBRARY_DIRS} )

else()
	find_path(NLOPT_INCLUDE_DIR nlopt.h
	          HINTS ${MUQ_NLOPT_DIR}/include
	          PATH_SUFFIXES nlopt NO_DEFAULT_PATH)

	find_library(NLOPT_LIBRARY NAMES nlopt nlopt_cxx
	             HINTS ${MUQ_NLOPT_DIR}/lib NO_DEFAULT_PATH)

	find_library(NLOPT_LIBRARY_STATIC NAMES ${library_prefix}nlopt.${static_library_suffix} ${library_prefix}nlopt_cxx.${static_library_suffix}
	             HINTS ${MUQ_NLOPT_DIR}/lib NO_DEFAULT_PATH)	 
endif()

set(NLOPT_LIBRARIES_STATIC ${NLOPT_LIBRARY_STATIC} )	
			 
set(NLOPT_LIBRARIES ${NLOPT_LIBRARY} )
set(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Nlopt  DEFAULT_MSG
                                  NLOPT_LIBRARY NLOPT_INCLUDE_DIR)

mark_as_advanced(NLOPT_INCLUDE_DIR NLOPT_LIBRARY )