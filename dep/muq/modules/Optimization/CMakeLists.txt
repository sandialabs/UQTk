add_subdirectory(src)
add_subdirectory(test)

if(MUQ_USE_PYTHON)
  add_subdirectory(python)
endif(MUQ_USE_PYTHON)
