project (UQTk)

SET(copy_FILES
  __init__.py
  pyuqtkarray_tools.py
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pyuqtkarray_tools.py
               ${CMAKE_CURRENT_BINARY_DIR}/pyuqtkarray_tools.py COPYONLY)

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/pyuqtkarray_tools)
