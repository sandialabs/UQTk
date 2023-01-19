project (UQTk)

SET(copy_FILES
 1_quadratic_one_data_set.ipynb
 2_quadratic_two_data_sets.ipynb
 3_quadratic_two_data_sets_disagree.ipynb
 4_quadratic_two_data_sets_different_nb_of_points.ipynb
 5_quadratic_with_optimal_beta.ipynb
 utils.py
 README
)

INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/dfi_app)