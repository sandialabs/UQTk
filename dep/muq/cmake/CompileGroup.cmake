# USAGE:
# CreateCompileGroup(
#     <group name>
#     <brief description>
#     <library_name>
#     <other compile group dependencies>
#     <required dependencies>
#     <optional dependencies>
#     <source file 1>
#     <source file 2>
#     ...
#     <source file N>
# )

# Specify a group of sources files, as well their dependencies, and ultimate target library.
function(CreateCompileGroup
         GROUP_NAME
	       DESCRIPTION
	       LIBRARY_NAME
	       GROUP_DEPENDENCIES
	       REQUIRED_DEPENDENCIES
	       OPTIONAL_DEPENDENCIES)

  option(MUQ_ENABLEGROUP_${GROUP_NAME} "Should the group ${GROUP_NAME} be compiled?" ${MUQ_ENABLEGROUP_DEFAULT})

  message(STATUS "MUQ_ENABLEGROUP_${GROUP_NAME} = ${MUQ_ENABLEGROUP_${GROUP_NAME}}")
  set(MUQ_GROUPS "${MUQ_GROUPS};${GROUP_NAME}" CACHE INTERNAL "A list of all of the muq groups.")

  set(${GROUP_NAME}_REQUIRES_GROUPS ${GROUP_DEPENDENCIES} CACHE INTERNAL "Other compile groups that are required by the ${GROUP_NAME} group.")

  set(${GROUP_NAME}_REQUIRES ${REQUIRED_DEPENDENCIES} CACHE INTERNAL "External packages required by the ${GROUP_NAME} compile group.")
  set(${GROUP_NAME}_DESIRES ${OPTIONAL_DEPENDENCIES} CACHE INTERNAL "External packages that the ${GROUP_NAME} compile group can exploit if they're available.")

  set(${GROUP_NAME}_LIBRARY ${LIBRARY_NAME} CACHE INTERNAL "The library this group will contribute to.")

  # Compute the path to the source file relative to the root directory
  string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/"  ""  RELATIVE_DIR ${CMAKE_CURRENT_LIST_DIR})

  set(SOURCES )
  foreach(source ${ARGN})
    list(APPEND SOURCES "${RELATIVE_DIR}/${source}")
  endforeach()

  set(${GROUP_NAME}_SOURCES ${SOURCES} CACHE INTERNAL "Source files included in the ${GROUP_NAME} compile group.")

endfunction(CreateCompileGroup)



# Used to add a collection of test files that are linked to a compile group specified with CreateCompileGroup
function(CreateTestGroup GROUP_NAME)

  set(MUQ_TEST_GROUPS "${MUQ_TEST_GROUPS};${GROUP_NAME}" CACHE INTERNAL "A list of all of the muq test groups.")

  # Compute the path to the source file relative to the root directory
  string(REPLACE "${CMAKE_SOURCE_DIR}/"  ""  RELATIVE_DIR ${CMAKE_CURRENT_LIST_DIR})

  set(SOURCES )
  foreach(source ${ARGN})
    list(APPEND SOURCES "${RELATIVE_DIR}/${source}")
  endforeach()

  set(${GROUP_NAME}_TEST_SOURCES ${SOURCES} CACHE INTERNAL "Source files with tests related to the ${GROUP_NAME} compile group.")

endfunction(CreateTestGroup)

# Used to add a collection of test files that are linked to a compile group specified with CreateCompileGroup
function(CreateParallelTestGroup GROUP_NAME)

  set(MUQ_PARALLEL_TEST_GROUPS "${MUQ_PARALLEL_TEST_GROUPS};${GROUP_NAME}" CACHE INTERNAL "A list of all of the muq parallel test groups.")

  # Compute the path to the source file relative to the root directory
  string(REPLACE "${CMAKE_SOURCE_DIR}/"  ""  RELATIVE_DIR ${CMAKE_CURRENT_LIST_DIR})

  set(SOURCES )
  foreach(source ${ARGN})
    list(APPEND SOURCES "${RELATIVE_DIR}/${source}")
  endforeach()

  set(${GROUP_NAME}_PARALLEL_TEST_SOURCES ${SOURCES} CACHE INTERNAL "Source files with parallel tests related to the ${GROUP_NAME} compile group.")

endfunction(CreateTestGroup)
