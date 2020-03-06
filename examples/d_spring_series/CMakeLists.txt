project (UQTk)

SET(copy_PY_FILES
    d_springs_tools.py
    run_d_springs.py
    README)

INSTALL(FILES ${copy_PY_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/d_spring_series)
