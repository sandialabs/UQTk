project (UQTk)

SET(copy_EX_FILES
    bare_bcs.ipynb
    readme.md)

INSTALL(FILES ${copy_EX_FILES}
        PERMISSIONS OWNER_WRITE OWNER_READ
        DESTINATION examples/bare_bcs)
