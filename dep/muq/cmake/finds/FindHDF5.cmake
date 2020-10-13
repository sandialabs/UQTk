
find_package(PkgConfig)

if(NOT DEFINED MUQ_HDF5_DIR)
	pkg_check_modules(PC_HDF5 QUIET libhdf5)
	set(HDF5_DEFINITIONS ${PC_HDF5_CFLAGS_OTHER})


	find_library(HDF5_LIBRARY NAMES hdf5
             HINTS ${PC_HDF5_LIBDIR} ${PC_HDF5_LIBRARY_DIRS} )

        if(HDF5_LIBRARY)
	    get_filename_component(TEMP_DIR ${HDF5_LIBRARY} DIRECTORY)
	    get_filename_component(PARENT_DIR ${TEMP_DIR} DIRECTORY)

	    find_library(HDF5_LIBRARY_STATIC NAMES ${library_prefix}hdf5.${static_library_suffix}
                 HINTS ${PARENT_DIR} ${PC_HDF5_LIBDIR} ${PC_HDF5_LIBRARY_DIRS} )

	    find_path(HDF5_INCLUDE_DIR NAMES hdf5.h
                 HINTS ${PARENT_DIR} ${PC_HDF5_INCLUDEDIR} ${PC_HDF5_INCLUDE_DIRS}
                 PATH_SUFFIXES lib hdf5 )


            pkg_check_modules(PC_HDF5HL QUIET libhdf5_hl)
            set(HDF5HL_DEFINITIONS ${PC_HDF5HL_CFLAGS_OTHER})

            find_library(HDF5HL_LIBRARY NAMES hdf5_hl
                         HINTS ${PARENT_DIR} ${PC_HDF5HL_LIBDIR} ${PC_HDF5HL_LIBRARY_DIRS}
                         PATH_SUFFIXES lib hdf5)

            find_library(HDF5HL_LIBRARY_STATIC NAMES ${library_prefix}hdf5_hl.${static_library_suffix}
                         HINTS ${PARENT_DIR} ${PC_HDF5HL_LIBDIR} ${PC_HDF5HL_LIBRARY_DIRS}
                         PATH_SUFFIXES lib hdf5)

        endif()
else()

	find_library(HDF5_LIBRARY NAMES hdf5
	             HINTS ${MUQ_HDF5_DIR}
                     PATH_SUFFIXES lib NO_DEFAULT_PATH)

        if(HDF5_LIBRARY)

	    get_filename_component(TEMP_DIR ${HDF5_LIBRARY} DIRECTORY)
	    get_filename_component(PARENT_DIR ${TEMP_DIR} DIRECTORY)

	    find_path(HDF5_INCLUDE_DIR NAMES hdf5.h
	              HINTS ${MUQ_HDF5_DIR} ${PARENT_DIR}
		      PATH_SUFFIXES include NO_DEFAULT_PATH)


            find_library(HDF5_LIBRARY_STATIC NAMES ${library_prefix}hdf5.${static_library_suffix}
                         HINTS ${MUQ_HDF5_DIR} ${PARENT_DIR}
                                 PATH_SUFFIXES lib NO_DEFAULT_PATH)

            find_library(HDF5HL_LIBRARY NAMES hdf5_hl
                         HINTS ${MUQ_HDF5_DIR}/lib ${PARENT_DIR}/lib NO_DEFAULT_PATH)

            find_library(HDF5HL_LIBRARY_STATIC NAMES ${library_prefix}hdf5_hl.${static_library_suffix}
                         HINTS ${MUQ_HDF5_DIR}/lib ${PARENT_DIR}/lib NO_DEFAULT_PATH)

        endif()

endif()

set(HDF5_LIBRARIES_STATIC ${HDF5_LIBRARY_STATIC} ${HDF5HL_LIBRARY_STATIC})

set(HDF5_LIBRARIES ${HDF5_LIBRARY} ${HDF5HL_LIBRARY})
set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HDF5  DEFAULT_MSG
                                  HDF5_LIBRARY HDF5_INCLUDE_DIR)

mark_as_advanced(HDF5_INCLUDE_DIR HDF5_LIBRARY )
