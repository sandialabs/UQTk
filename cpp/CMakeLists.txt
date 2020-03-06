project (UQTk)

link_directories(${CMAKE_SUNDIALS_DIR}/lib)
if( BUILD_SUNDIALS)
  link_directories("${PROJECT_BINARY_DIR}/../dep/sundials/src/cvode")
  link_directories("${PROJECT_BINARY_DIR}/../dep/sundials/src/nvector/serial")
  link_directories("${PROJECT_BINARY_DIR}/../dep/sundials/src/sunlinsol/dense")
  link_directories("${PROJECT_BINARY_DIR}/../dep/sundials/src/sunmatrix/dense")
endif()

add_subdirectory (lib)
add_subdirectory (app)
add_subdirectory (tests)
