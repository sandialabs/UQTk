# Using Pybind to build the python module

find_package(pybind11 REQUIRED)
pybind11_add_module(pyuqtkarray pyuqtkarray.cpp)

INSTALL(TARGETS pyuqtkarray DESTINATION PyUQTk/)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${outdir}PyUQTk/pyuqtkarray/pyuqtkarray.py
	      DESTINATION PyUQTk)
