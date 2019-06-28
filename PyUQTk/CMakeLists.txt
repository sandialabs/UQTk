project (UQTk) 

add_subdirectory (inference)
add_subdirectory (plotting)
add_subdirectory (sens)
add_subdirectory (multirun)
add_subdirectory (utils)

SET(copy_FILES
  __init__.py
  )
INSTALL(FILES ${copy_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION PyUQTk/
)

if ("${PyUQTk}" STREQUAL "ON")
  include(numpy.cmake)

  add_subdirectory (uqtkarray)
  add_subdirectory (quad)
  add_subdirectory (tools)
  #add_subdirectory (kle)
  add_subdirectory (pce)
  add_subdirectory (bcs)
  # add_subdirectory (mcmc)
  #add_subdirectory (dfi)

  add_subdirectory(pytests)
endif()
