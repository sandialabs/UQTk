project (UQTk)

add_subdirectory (line)

SET(copy_DOC_FILES
    README)

INSTALL(FILES ${copy_DOC_FILES}
        PERMISSIONS OWNER_WRITE OWNER_READ
        DESTINATION examples/dfi)
