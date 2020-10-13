cmake_minimum_required (VERSION 2.8.12)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/finds)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/checks)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/builds)

INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckCXXSourceRuns)
INCLUDE(CheckFunctionExists)

# define the project name
project(MUQ)
set(${CMAKE_PROJECT_NAME}_MAJOR_VERSION 0)
set(${CMAKE_PROJECT_NAME}_MINOR_VERSION 2)
set(${CMAKE_PROJECT_NAME}_PATCH_VERSION 0)
set(${CMAKE_PROJECT_NAME}_VERSION ${${CMAKE_PROJECT_NAME}_MAJOR_VERSION}.${${CMAKE_PROJECT_NAME}_MINOR_VERSION}.${${CMAKE_PROJECT_NAME}_PATCH_VERSION})

message(STATUS "\n\nThis is MUQ, the MIT UQ Library!\n\n")

# Make the external project installation directory if it doesn't exist already
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/include)
file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/muq_external/lib)

###############################################
# Compiler configuration
###############################################
INCLUDE(SetupCompiler)
INCLUDE(CompilerCheck)
include_directories(.)

###############################################
# User options
###############################################
INCLUDE(OptionsSetup)

###############################################
# Find source files and calculate required dependencies
###############################################
INCLUDE(CompileGroup)

# The ProcesCompileGroups function will add modules/CMakeLists.txt and then process any computed Compile/Test groups
INCLUDE(ProcessCompileGroups)

###############################################
# Manage dependencies
###############################################
INCLUDE(DependencySearch)

include_directories(. ./external/include)

# write a configure file containing information about what libraries are available
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/MUQ/config.h)

###############################################
# set module dependencies
###############################################
INCLUDE(DependencyCheck)

###############################################
# Create BUILD targets for cmake
###############################################

# Remove any compile groups with unsatisfied dependencies
INCLUDE(UpdateCompileGroups)

# Set up libraries and test executables based on the compile groups
INCLUDE(SetupCompileTargets)
INCLUDE(SetupTestTargets)

###############################################################
# INSTALLATION INFORMATION
###############################################################
get_directory_property(MUQ_COMPILE_DEFINITIONS COMPILE_DEFINITIONS)
message(STATUS "Storing MUQ_COMPILE_DEFINITIONS = ${MUQ_COMPILE_DEFINITIONS}")
message(STATUS "MUQ_LIBRARIES = ${MUQ_LIBRARIES}")

set(MUQ_LIBRARIES_CACHE ${MUQ_LIBRARIES} CACHE INTERNAL "List of MUQ libraries for linking in parent projects.")
set(MUQ_LINK_LIBRARIES_CACHE ${MUQ_LINK_LIBRARIES} CACHE INTERNAL "List of external libraries needed by MUQ for use in parent projects.")

set(MUQ_EXPORT_LIBRARIES "")
FOREACH(BUILT_LIB ${MUQ_LIBRARIES})
	list(APPEND MUQ_EXPORT_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}${BUILT_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
ENDFOREACH()
message(STATUS "MUQ_EXPORT_LIBRARIES = ${MUQ_EXPORT_LIBRARIES}")

#export(PACKAGE ${CMAKE_PROJECT_NAME})

# Create the MUQConfig.cmake file
configure_file(cmake/MUQConfig.cmake.in "${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)
configure_file(cmake/MUQConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${CMAKE_PROJECT_NAME}Config.cmake" @ONLY)

# Install the MUQConfig.cmake file
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${CMAKE_PROJECT_NAME}Config.cmake"
        DESTINATION "${CMAKE_INSTALL_PREFIX}/CMake/${CMAKE_PROJECT_NAME}"
	    COMPONENT dev)

# Install the pymuq libraries
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/Modeling/python/pymuqModeling/"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pymuqModeling"
	    COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT ${CMAKE_PROJECT_NAME}Depends
        DESTINATION "${CMAKE_INSTALL_PREFIX}/CMake/${CMAKE_PROJECT_NAME}"
	    COMPONENT dev)

# install all the MUQ headers -- may be better off installing for each package individually
install(DIRECTORY MUQ
        DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
	      FILES_MATCHING PATTERN "*.h")

# install the python libraries -- Is there a way to avoid copying this for each library?
#install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/Approximation/python/pymuqApproximation/"
#        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pymuqApproximation"
#        PATTERN "wrappers/" EXCLUDE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${INSTALL_LIB_DIR})

###############################################
# Write summary to terminal
###############################################

# write a summary of the build configuration
include(MakeSummaryFile)

FILE(READ ${CMAKE_BINARY_DIR}/summary.log MUQ_LOG_SUMMARY)
MESSAGE("${MUQ_LOG_SUMMARY}")


###############################################
# Add a target to build documentation
###############################################
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	message(STATUS "To build documentation, run 'make doc'")

	option(MUQDOC_INCLUDE_ANALYTICS OFF)
	if(MUQDOC_INCLUDE_ANALYTICS)
		message(STATUS "Will include google analytics header information in doxygen output.")
		set(ANALYTICS_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/documentation/doxFiles/analytics_header.html)
	endif()

  message(STATUS "Extracting parameter values from cpp files.  Results are stored in ${CMAKE_BINARY_DIR}/ExtractedParameters.dox")
  message(STATUS "python ${CMAKE_CURRENT_SOURCE_DIR}/documentation/doxFiles/FindPtreeParameters.py ${CMAKE_CURRENT_SOURCE_DIR}")
  execute_process(COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/documentation/doxFiles/FindPtreeParameters.py ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_FILE "${CMAKE_BINARY_DIR}/ParameterDescriptions.dox")

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/documentation/muq.doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/muq.doxyfile @ONLY)
    #configure_file(${CMAKE_CURRENT_SOURCE_DIR}/documentation/ParameterDescriptions.dox.in ${CMAKE_CURRENT_BINARY_DIR}/ParameterDescriptions.dox @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/muq.doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )


endif(DOXYGEN_FOUND)
