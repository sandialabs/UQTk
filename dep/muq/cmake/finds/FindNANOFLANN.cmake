find_package(PkgConfig)
if(NOT DEFINED MUQ_NANOFLANN_DIR)
	
	pkg_check_modules(PC_NANOFLANN QUIET NANOFLANN)
	set(NANOFLANN_DEFINITIONS ${PC_NANOFLANN_CFLAGS_OTHER})

	find_path(NANOFLANN_INCLUDE_DIR 
    		  HINTS ${PC_NANOFLANN_INCLUDEDIR} ${PC_NANOFLANN_INCLUDE_DIRS}
          	  PATH_SUFFIXES nanoflann)

else()

	find_path(NANOFLANN_INCLUDE_DIR 
	          HINTS ${MUQ_NANOFLANN_DIR}
			  PATH_SUFFIXES nanoflann NO_DEFAULT_PATH)
					
endif()
set(NANOFLANN_INCLUDE_DIRS ${NANOFLANN_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(nanoflann DEFAULT_MSG NANOFLANN_INCLUDE_DIR)

mark_as_advanced(NANOFLANN_INCLUDE_DIR)
