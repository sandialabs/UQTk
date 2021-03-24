project (UQTk)

SET(copy_FILES
  __init__.py
  pce_tools.py
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pce_tools.py
               ${CMAKE_CURRENT_BINARY_DIR}/pce_tools.py COPYONLY)

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/pce_tools)
