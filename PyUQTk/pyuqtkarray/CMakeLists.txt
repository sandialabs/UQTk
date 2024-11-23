# Using Pybind to build the Python module
find_package(pybind11 REQUIRED)

# Include directories for UQTk headers
include_directories(../../cpp/lib/array) # Array directory
include_directories(../../cpp/lib/include) # Include directory

# Add the Pybind11 module
pybind11_add_module(_uqtkarray pyuqtkarray.cpp)

# Link required libraries
target_link_libraries(_uqtkarray PUBLIC uqtk blas lapack)

# Files to be copied
set(copy_FILES
  pyuqtkarray_tools.py
)

# Install Python files
install(FILES ${copy_FILES} DESTINATION PyUQTk/)

# Ensure Python tools script is available in the build directory
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/pyuqtkarray_tools.py
  ${CMAKE_CURRENT_BINARY_DIR}/pyuqtkarray_tools.py COPYONLY
)

# Install the compiled Python module (_uqtkarray)
install(TARGETS _uqtkarray DESTINATION PyUQTk/)