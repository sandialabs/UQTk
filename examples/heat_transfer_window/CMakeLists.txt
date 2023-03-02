project (UQTk)

SET(copy_PY_FILES
    heat_transfer_window.ipynb
    window_tools.py
    img_window_1.png
    img_window_2.png	
    README)

INSTALL(FILES ${copy_PY_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/heat_transfer_window)
