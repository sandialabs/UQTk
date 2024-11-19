# Use modern Python detection (requires CMake 3.12+)
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Check if NumPy was found
if (Python3_NumPy_FOUND)
    # Extract NumPy include directory and version
    set(NUMPY_INCLUDE_DIR ${Python3_NumPy_INCLUDE_DIRS})
    set(NUMPY_VERSION ${Python3_NumPy_VERSION})

    # Debugging output
    message(STATUS "Found NumPy:")
    message(STATUS "  Include Directory: ${NUMPY_INCLUDE_DIR}")
    message(STATUS "  Version: ${NUMPY_VERSION}")
else()
    # Handle error if NumPy is not found
    message(FATAL_ERROR "NumPy not found. Please ensure NumPy is installed in your Python3 environment.")
endif()

# Include NumPy headers
include_directories(${NUMPY_INCLUDE_DIR})