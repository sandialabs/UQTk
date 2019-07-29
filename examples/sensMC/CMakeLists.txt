project (UQTk)

SET(copy_SCRIPT_FILES
  run.sh
  model.sh
  genRndSpls.sh
  )

INSTALL(
        FILES ${copy_SCRIPT_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/sensMC
)

