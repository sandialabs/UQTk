project (UQTk)

SET(copy_FILES
  muqPCE2uqtk.cpp
  genzGP.py
)

# muqPCE2uqtk.x
link_directories(${CMAKE_SUNDIALS_DIR}/lib)
add_executable (muqPCE2uqtk.x muqPCE2uqtk.cpp)

# UQTk libraries
target_link_libraries (muqPCE2uqtk.x uqtkpce  )
target_link_libraries (muqPCE2uqtk.x uqtkmcmc )
target_link_libraries (muqPCE2uqtk.x uqtkquad )
target_link_libraries (muqPCE2uqtk.x uqtktools)
target_link_libraries (muqPCE2uqtk.x uqtkarray )
target_link_libraries (muqPCE2uqtk.x uqtkxmlutils )

target_link_libraries (muqPCE2uqtk.x depdsfmt )
target_link_libraries (muqPCE2uqtk.x deplbfgs )
target_link_libraries (muqPCE2uqtk.x sundials_cvode)
target_link_libraries (muqPCE2uqtk.x sundials_nvecserial)
target_link_libraries (muqPCE2uqtk.x sundials_sunlinsoldense)
target_link_libraries (muqPCE2uqtk.x sundials_sunmatrixdense)
target_link_libraries (muqPCE2uqtk.x depslatec )
target_link_libraries (muqPCE2uqtk.x deplapack )
target_link_libraries (muqPCE2uqtk.x depblas )
target_link_libraries (muqPCE2uqtk.x depfigtree )
target_link_libraries (muqPCE2uqtk.x depann   )
target_link_libraries (muqPCE2uqtk.x expat )

# MUQ libraries
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libhdf5_hl.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libhdf5.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/lib/libmuqModelling.dylib)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/lib/libmuqApproximation.dylib)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/lib/libmuqUtilities.dylib)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libflann_cpp_s.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libsundials_cvodes.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libsundials_idas.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libsundials_kinsol.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libsundials_nvecserial.a)
#target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libhdf5_hl.a)
#target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libhdf5.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libboost_filesystem.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libboost_system.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libboost_graph.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libboost_serialization.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libboost_regex.a)
target_link_libraries (muqPCE2uqtk.x ${PATH2MUQ}/muq_external/lib/libboost_wserialization.a)

# Link fortran libraries
# using GCC
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-maybe-uninitialized -Wno-sign-compare")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-local-typedefs")
if ("${GnuLibPath}" STREQUAL "")
  target_link_libraries (muqPCE2uqtk.x gfortran stdc++)
else()
  target_link_libraries (muqPCE2uqtk.x  ${GnuLibPath}/libgfortran.a ${GnuLibPath}/libquadmath.a stdc++)
endif()

include_directories(../../examples/line_infer)
include_directories(../../cpp/lib/include)
include_directories(../../cpp/lib/array)
include_directories(../../cpp/lib/pce)
include_directories(../../cpp/lib/quad)
include_directories(../../cpp/lib/mcmc)
include_directories(../../cpp/lib/tools)
include_directories(../../cpp/lib/xmlutils)

include_directories(../../dep/dsfmt)
include_directories(../../dep/figtree)
include_directories(../../dep/cvode-2.7.0/include)
include_directories("${PROJECT_BINARY_DIR}/../../dep/cvode-2.7.0/include")

include_directories(${PATH2MUQ}/include)
include_directories(${PATH2MUQ}/muq_external/include/eigen3)
include_directories(${PATH2MUQ}/muq_external/include)
include_directories(/usr/local/include/eigen3)
include_directories(/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7)

INSTALL(TARGETS muqPCE2uqtk.x DESTINATION examples/muq)
INSTALL(FILES ${copy_FILES} DESTINATION examples/muq)
