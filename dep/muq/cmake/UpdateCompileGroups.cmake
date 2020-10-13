
# Turn off any compile groups with unfulfilled external library dependencies
foreach(group ${MUQ_GROUPS})
    foreach(depend ${${group}_REQUIRES})
        #if(NOT MUQ_USE_${depend})
        #    message(STATUS "Turning off ${group} because of a missing library dependency.")
        #    set(MUQ_ENABLEGROUP_${group} OFF)
        #endif()
    endforeach()
endforeach()

# Turn off any compile groups whose upstream compile groups are off
foreach(group ${MUQ_GROUPS})
    foreach(depend ${${group}_REQUIRES_GROUPS})
        if(MUQ_ENABLEGROUP_${group} AND NOT MUQ_ENABLEGROUP_${depend})
            message(STATUS "Turning off ${group} because of a missing upstream group dependency (${depend}).")
            set(MUQ_ENABLEGROUP_${group} OFF)
        endif()
    endforeach()
endforeach()

# Create a list of all MUQ libraries to build
set(MUQ_TARGETS )
foreach(group ${MUQ_GROUPS})
    if(MUQ_ENABLEGROUP_${group})
        list(APPEND MUQ_TARGETS ${${group}_LIBRARY})
    endif()
endforeach()
list(REMOVE_DUPLICATES MUQ_TARGETS)

# Set up the source for each target library
foreach(target ${MUQ_TARGETS})
    set(${target}_SOURCES )

    foreach(group ${MUQ_GROUPS})
        if(MUQ_ENABLEGROUP_${group})
	    if(${${group}_LIBRARY} STREQUAL ${target})

                # Check to see if a group has any source (e.g., *.cpp) files.  Flag it as something that will be built if it does.
	        list(LENGTH ${group}_SOURCES sources_length)
		if(sources_length GREATER 0)
	            set(${group}_IS_COMPILED ON CACHE INTERNAL "Whether or not the ${group} is used in any library.")
		endif()

	        list(APPEND ${target}_SOURCES ${${group}_SOURCES})
	    endif()
	endif()
    endforeach()

    if(${target}_SOURCES)
        list(REMOVE_DUPLICATES ${target}_SOURCES)
    endif()

endforeach()
