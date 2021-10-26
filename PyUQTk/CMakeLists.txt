project (UQTk)

add_subdirectory (inference)
add_subdirectory (plotting)
add_subdirectory (sens)
add_subdirectory (multirun)
add_subdirectory (utils)
if ("${PyUQTk}" STREQUAL "ON")
  add_subdirectory (PyPCE)
  add_subdirectory (pyuqtkarray_tools)
  add_subdirectory (bcs_ext)
endif()

SET(copy_FILES
  __init__.py
  uqtkarray.py
  tools.py
  quad.py
  pce.py
  bcs.py
  kle.py
  )
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

  add_subdirectory(pytests)
endif()

# Copy over README file too
INSTALL(FILES README
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION PyUQTk/)
