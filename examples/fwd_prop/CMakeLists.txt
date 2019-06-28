project (UQTk)

SET(copy_PY_FILES
    rad_heat_transfer_atm_pce.py
    heat_transfer_pce_tools.py)

INSTALL(FILES ${copy_PY_FILES} 
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ	
        DESTINATION examples/fwd_prop)
