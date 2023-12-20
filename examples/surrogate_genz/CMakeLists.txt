project (UQTk)

SET(copy_PY_FILES
    surrogate_genz-BCS.ipynb
    surrogate_genz-Galerkin.ipynb
    surrogate_genz-Regression.ipynb
    surrogate_genz.ipynb	
    README.md)

INSTALL(FILES ${copy_PY_FILES}
        PERMISSIONS OWNER_WRITE OWNER_READ
        DESTINATION examples/surrogate_genz)
