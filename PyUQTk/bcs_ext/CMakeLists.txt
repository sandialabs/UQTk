project (UQTk)

include_directories(../pyuqtkarray)
include_directories(../pce)
include_directories(../tools)
include_directories(../pyuqtkarray_tools)
include_directories(../bcs)

SET(copy_FILES
  __init__.py
  bcs_ext.py
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/bcs_ext.py
               ${CMAKE_CURRENT_BINARY_DIR}/bcs_ext.py COPYONLY)

INSTALL(FILES ${copy_FILES} DESTINATION PyUQTk/bcs_ext)