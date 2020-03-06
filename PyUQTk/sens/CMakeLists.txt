project (UQTk)

SET(copy_FILES
  __init__.py
  gsalib.py
  )

INSTALL(FILES ${copy_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION PyUQTk/sens)
