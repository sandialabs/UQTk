project (UQTk)

SET(copy_FILES
  evidence.py
  fit.py
  get_data.py
  post.py
  tools.py
  graph_tools.py)


INSTALL(FILES ${copy_FILES}
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
        DESTINATION examples/polynomial/scripts)
