project (UQTk)

SET(copy_FILES
 uq_pc.py
 model.py
 plot.py
 plot_prep.py
 join_results.py
 prepare_inpc.py
 generate_inputsamples.py
 get_modelpc.py
 model_sens.x
 transpose_file.x
 scale.x
 col_ave.x
 getrange.x
 example_0.x
 example_1.x
 example_2.x
 example_3.x
 )


INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/uqpc)
