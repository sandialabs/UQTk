add_subdirectory (gen_mi)
add_subdirectory (generate_quad)
add_subdirectory (model_inf)
add_subdirectory (pce_eval)
add_subdirectory (pce_quad)
add_subdirectory (pce_resp)
add_subdirectory (pce_rv)
add_subdirectory (pce_sens)
add_subdirectory (pdf_cl)
add_subdirectory (sens)
add_subdirectory (regression)
add_subdirectory (gp_regr)
add_subdirectory (gkpSparse)
add_subdirectory (lr_regr)
add_subdirectory (lr_eval)

INSTALL(FILES README
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION bin/)
