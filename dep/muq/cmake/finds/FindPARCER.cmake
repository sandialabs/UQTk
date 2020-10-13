
find_package(PkgConfig)

if(NOT DEFINED MUQ_PARCER_DIR)

	pkg_check_modules(PC_PARCER QUIET libparcer)
	set(PARCER_DEFINITIONS ${PC_PARCER_CFLAGS_OTHER})

	find_path(PARCER_INCLUDE_DIR parcer/Communicator.h
					HINTS ${PC_PARCER_INCLUDEDIR} ${PC_PARCER_INCLUDE_DIRS}
					PATH_SUFFIXES parcer )

	find_library(PARCER_LIBRARY NAMES parcer
             HINTS ${PC_PARCER_LIBDIR} ${PC_PARCER_LIBRARY_DIRS} )

	find_library(PARCER_LIBRARY_STATIC NAMES ${library_prefix}parcer.${static_library_suffix}
             HINTS ${PC_PARCER_LIBDIR} ${PC_PARCER_LIBRARY_DIRS} )

else()
	find_path(PARCER_INCLUDE_DIR parcer/Communicator.h
	          HINTS ${MUQ_PARCER_DIR}/include
	          PATH_SUFFIXES parcer NO_DEFAULT_PATH)

	find_library(PARCER_LIBRARY NAMES parcer
	             HINTS ${MUQ_PARCER_DIR}/lib NO_DEFAULT_PATH)

endif()

set(PARCER_LIBRARIES ${PARCER_LIBRARY} )
set(PARCER_INCLUDE_DIRS ${PARCER_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(ParCer  DEFAULT_MSG
                                  PARCER_LIBRARY PARCER_INCLUDE_DIR)

mark_as_advanced(PARCER_INCLUDE_DIR PARCER_LIBRARY )
