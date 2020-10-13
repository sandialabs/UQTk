cmake_minimum_required (VERSION 2.8)

find_package(MUQ_PRIVATE REQUIRED)

include_directories(${MUQ_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})
set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})


add_executable(MonotoneRegression MonotoneRegression.cpp)

target_link_libraries(MonotoneRegression ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
