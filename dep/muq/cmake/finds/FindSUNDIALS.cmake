

if(NOT DEFINED MUQ_SUNDIALS_DIR)

  find_path(SUNDIALS_INCLUDE_DIR cvodes/cvodes.h
            HINTS ${SUNDIALS_INCLUDE_DIRS} /usr/local/include/ /usr/include/)
  find_library(CVODES_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_cvodes${CMAKE_SHARED_LIBRARY_SUFFIX}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(IDAS_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_idas${CMAKE_SHARED_LIBRARY_SUFFIX}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(KINSOL_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_kinsol${CMAKE_SHARED_LIBRARY_SUFFIX}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(NVEC_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_nvecserial${CMAKE_SHARED_LIBRARY_SUFFIX}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(NVEC_PARALLEL_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_nvecparallel${CMAKE_SHARED_LIBRARY_SUFFIX}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)

  find_library(CVODES_LIBRARY_STATIC NAMES ${library_prefix}sundials_cvodes${static_library_suffix}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(IDAS_LIBRARY_STATIC NAMES ${library_prefix}sundials_idas${static_library_suffix}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(KINSOL_LIBRARY_STATIC NAMES ${library_prefix}sundials_kinsol${static_library_suffix}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(NVEC_LIBRARY_STATIC NAMES ${library_prefix}sundials_nvecserial${static_library_suffix}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
  find_library(NVEC_PARALLEL_LIBRARY_STATIC NAMES ${library_prefix}sundials_nvecparallel${static_library_suffix}
               HINTS ${SUNDIALS_LIBRARY_DIRS} /usr/local/lib/ /usr/lib/)
else()

	find_path(SUNDIALS_INCLUDE_DIR cvodes/cvodes.h
	          HINTS ${MUQ_SUNDIALS_DIR}/include/ NO_DEFAULT_PATH)
	find_library(CVODES_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_cvodes${CMAKE_SHARED_LIBRARY_SUFFIX}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
	find_library(IDAS_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_idas${CMAKE_SHARED_LIBRARY_SUFFIX}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
	find_library(KINSOL_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_kinsol${CMAKE_SHARED_LIBRARY_SUFFIX}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
	find_library(NVEC_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_nvecserial${CMAKE_SHARED_LIBRARY_SUFFIX}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
  find_library(NVEC_PARALLEL_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}sundials_nvecparallel${CMAKE_SHARED_LIBRARY_SUFFIX}
             	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)

	find_library(CVODES_LIBRARY_STATIC NAMES ${library_prefix}sundials_cvodes${static_library_suffix}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
	find_library(IDAS_LIBRARY_STATIC NAMES ${library_prefix}sundials_idas${static_library_suffix}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
	find_library(KINSOL_LIBRARY_STATIC NAMES ${library_prefix}sundials_kinsol${static_library_suffix}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
	find_library(NVEC_LIBRARY_STATIC NAMES ${library_prefix}sundials_nvecserial${static_library_suffix}
	             HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
  find_library(NVEC_PARALLEL_LIBRARY_STATIC NAMES ${library_prefix}sundials_nvecparallel${static_library_suffix}
                            HINTS ${MUQ_SUNDIALS_DIR}/lib/ NO_DEFAULT_PATH)
endif()

set(SUNDIALS_LIBRARY ${CVODES_LIBRARY} ${IDAS_LIBRARY} ${KINSOL_LIBRARY} ${NVEC_LIBRARY} ${NVEC_PARALLEL_LIBRARY})
set(SUNDIALS_LIBRARIES ${SUNDIALS_LIBRARY})

set(SUNDIALS_LIBRARY_STATIC ${CVODES_LIBRARY_STATIC} ${IDAS_LIBRARY_STATIC} ${KINSOL_LIBRARY_STATIC} ${NVEC_LIBRARY_STATIC} ${NVEC_PARALLEL_LIBRARY_STATIC})
set(SUNDIALS_LIBRARIES_STATIC ${SUNDIALS_LIBRARY_STATIC})

set(SUNDIALS_INCLUDE_DIRS ${SUNDIALS_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS  DEFAULT_MSG
	                              SUNDIALS_LIBRARY SUNDIALS_INCLUDE_DIR)
mark_as_advanced(SUNDIALS_INCLUDE_DIR SUNDIALS_LIBRARY )

if( CVODES_LIBRARY OR CVODES_LIBRARY_STATIC )
    set(MUQ_HAS_SUNDIALS 1)
else()
    set(MUQ_HAS_SUNDIALS 0)
endif()
