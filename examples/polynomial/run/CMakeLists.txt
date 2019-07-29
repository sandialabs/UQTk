project (UQTk)

SET(copy_FILES
  full_run.sh
  input.xml)


INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/polynomial/run)
