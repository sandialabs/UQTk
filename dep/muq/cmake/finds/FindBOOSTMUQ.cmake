
find_package(PkgConfig)

set(BOOST_MIN_VERSION "1.56.0")

if(NOT DEFINED MUQ_BOOST_DIR)

	unset(Boost_LIBRARIES)
	unset(Boost_INCLUDE_DIR)
	unset(Boost_LIBRARY_DIRS)
	set(Boost_USE_STATIC_LIBS ON)

	find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem graph)

	IF(Boost_FOUND)
		set(BOOST_LIBRARIES_STATIC ${Boost_LIBRARIES})
	endif()

	unset(Boost_LIBRARIES)
	unset(Boost_INCLUDE_DIR)
	unset(Boost_LIBRARY_DIRS)
	unset(Boost_USE_STATIC_LIBS)

	find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS system filesystem graph)

	IF(Boost_FOUND)
		set(BOOST_LIBRARY ${Boost_LIBRARIES})
		set(BOOST_INCLUDE_DIR ${Boost_INCLUDE_DIR})
	endif()

else()

	find_path(BOOST_INCLUDE_DIR boost/property_tree/ptree.hpp
	          HINTS ${MUQ_BOOST_DIR}/include ${MUQ_BOOST_DIR}
	          PATH_SUFFIXES boost NO_DEFAULT_PATH)

	find_library(BOOST_SYSTEM_LIBRARY_STATIC NAMES ${library_prefix}boost_system.${static_library_suffix}
	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
	find_library(BOOST_FILESYSTEM_LIBRARY_STATIC NAMES ${library_prefix}boost_filesystem.${static_library_suffix}
	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
 	find_library(BOOST_GRAPH_LIBRARY_STATIC NAMES ${library_prefix}boost_graph.${static_library_suffix}
 	             HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)

 	find_library(BOOST_SYSTEM_LIBRARY NAMES boost_system
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
 	find_library(BOOST_FILESYSTEM_LIBRARY NAMES boost_filesystem
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)
 	find_library(BOOST_GRAPH_LIBRARY NAMES boost_graph
 		 	     HINTS ${MUQ_BOOST_DIR}/lib ${MUQ_BOOST_DIR}/stage/lib NO_DEFAULT_PATH)

	set(BOOST_LIBRARY ${BOOST_SYSTEM_LIBRARY} ${BOOST_FILESYSTEM_LIBRARY} ${BOOST_GRAPH_LIBRARY})
	set(BOOST_LIBRARY_STATIC ${BOOST_SYSTEM_LIBRARY_STATIC} ${BOOST_FILESYSTEM_LIBRARY_STATIC} ${BOOST_GRAPH_LIBRARY_STATIC})
endif()


set(BOOST_INCLUDE_DIRS ${BOOST_INCLUDE_DIR})
set(BOOST_LIBRARIES ${BOOST_LIBRARY})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(BOOST  DEFAULT_MSG
                                  BOOST_LIBRARY BOOST_INCLUDE_DIR)

mark_as_advanced(BOOST_INCLUDE_DIR BOOST_LIBRARY)
