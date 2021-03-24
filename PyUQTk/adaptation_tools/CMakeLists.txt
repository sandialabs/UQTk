project (UQTk)

SET(copy_FILES
  __init__.py
  adaptation_tools.py
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/adaptation_tools.py
               ${CMAKE_CURRENT_BINARY_DIR}/adaptation_tools.py COPYONLY)

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/adaptation_tools)
