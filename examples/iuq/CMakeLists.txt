project (UQTk)

SET(copy_FILES
 workflow.x
 prep_model1.x
 prep_data1.x
 prep_model2.x
 prep_data2.x
 prep_calib.x
 run_infer.x
 postp_infer.x
 get_postpred.py
 get_postsam.py
 plot_sens.py
 plot_pdfs.py
 plot_prpost.py
 plot_fit1d.py
 plot_shade.py
 README
 )


INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/iuq)
