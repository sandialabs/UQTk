project (UQTk)

SET(copy_EXE_FILES
    spl_qmc.x
   )

# spl_qmc.x
add_executable (spl_qmc.x spl_qmc.cpp)
target_link_libraries (spl_qmc.x uqtksampling uqtkarray)

include_directories(../../cpp/lib/include )
include_directories(../../cpp/lib/sampling)
include_directories(../../cpp/lib/array   )
include_directories(../../cpp/lib/tools   )
include_directories(../../dep/dsfmt       )
include_directories(../../dep/figtree     )

INSTALL(TARGETS ${copy_EXE_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/sampling_qmc
)
