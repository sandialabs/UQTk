project (UQTk)

SET(copy_PY_FILES
    5D_fwd_prop.py
    fwd_prop_tools.py
    highd_sparse.py
    tools_conductance_dp_pce_wrad.py)

INSTALL(FILES ${copy_PY_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ	
        DESTINATION examples/window)
