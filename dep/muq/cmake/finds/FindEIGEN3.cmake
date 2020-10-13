find_package(PkgConfig)
if(NOT DEFINED MUQ_EIGEN3_DIR)

	pkg_check_modules(PC_EIGEN3 QUIET EIGEN3)
	set(EIGEN3_DEFINITIONS ${PC_EIGEN3_CFLAGS_OTHER})

	find_path(EIGEN3_INCLUDE_DIR Eigen/Core
    		  HINTS ${PC_EIGEN3_INCLUDEDIR} ${PC_EIGEN3_INCLUDE_DIRS}
          	  PATH_SUFFIXES eigen3)

else()

	find_path(EIGEN3_INCLUDE_DIR Eigen/Core
	          HINTS ${MUQ_EIGEN3_DIR}
			  PATH_SUFFIXES eigen3 NO_DEFAULT_PATH)

endif()
set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIR)

mark_as_advanced(EIGEN3_INCLUDE_DIR)
