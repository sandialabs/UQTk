project (UQTk)

add_subdirectory (scripts)
add_subdirectory (run)

# Copy over README file too
INSTALL(FILES README
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/polynomial)
