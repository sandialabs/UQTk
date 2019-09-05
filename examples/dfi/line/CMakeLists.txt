project (UQTk)

SET(copy_FILES
  main.cpp
  Makefile
  userFunctions_line.cpp
  )

INSTALL(
        FILES ${copy_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/dfi/line
)

