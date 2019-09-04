project (UQTk)

SET(copy_FILES
  main.cpp
  userFunctions_line.cpp
  Makefile
  )

INSTALL(FILES ${copy_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
        DESTINATION examples/dfi_line)
