project (UQTk)

add_subdirectory (inference)
add_subdirectory (plotting)
add_subdirectory (sens)
add_subdirectory (multirun)
add_subdirectory (utils)

if ("${PyUQTk}" STREQUAL "ON")
  add_subdirectory (PyPCE)
endif()

SET(copy_FILES
  __init__.py
  uqtkarray.py
  tools.py
  quad.py
  pce.py
  bcs.py
  kle.py
  lowrank.py
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/uqtkarray.py
               ${CMAKE_CURRENT_BINARY_DIR}/uqtkarray.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tools.py
               ${CMAKE_CURRENT_BINARY_DIR}/tools.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/quad.py
               ${CMAKE_CURRENT_BINARY_DIR}/quad.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pce.py
               ${CMAKE_CURRENT_BINARY_DIR}/pce.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bcs.py
               ${CMAKE_CURRENT_BINARY_DIR}/bcs.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/kle.py
               ${CMAKE_CURRENT_BINARY_DIR}/kle.py COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lowrank.py
               ${CMAKE_CURRENT_BINARY_DIR}/lowrank.py COPYONLY)

INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION PyUQTk/
)

if ("${PyUQTk}" STREQUAL "ON")
  include(numpy.cmake)

  add_subdirectory (pyuqtkarray)
  add_subdirectory (quad)
  add_subdirectory (tools)
  add_subdirectory (kle)
  add_subdirectory (pce)
  add_subdirectory (bcs)
  add_subdirectory (lowrank)

  add_subdirectory(pytests)
endif()

# Copy over README file too
INSTALL(FILES README
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION PyUQTk/)
