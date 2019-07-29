project (UQTk)

SET(copy_PY_FILES
    full_quad.py
    sparse_quad.py
    quad_tools.py
)

INSTALL(
        FILES ${copy_PY_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ	
        DESTINATION examples/num_integ
)
