find_package(PkgConfig)
if(NOT DEFINED MUQ_STANMATH_DIR)

	pkg_check_modules(PC_EIGEN3 QUIET STANMATH)
	set(EIGEN3_DEFINITIONS ${PC_STANMATH_CFLAGS_OTHER})

	find_path(EIGEN3_INCLUDE_DIR Eigen/Core
    		  HINTS ${PC_STANMATH_INCLUDEDIR} ${PC_STANMATH_INCLUDE_DIRS}
          	  PATH_SUFFIXES eigen3)

else()

	find_path(STANMATH_INCLUDE_DIR stan/math.hpp
	          HINTS ${MUQ_STANMATH_DIR}
			  PATH_SUFFIXES stan NO_DEFAULT_PATH)

endif()
set(STANMATH_INCLUDE_DIRS ${STANMATH_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(STANMATH DEFAULT_MSG STANMATH_INCLUDE_DIR)

mark_as_advanced(STANMATH_INCLUDE_DIR)
