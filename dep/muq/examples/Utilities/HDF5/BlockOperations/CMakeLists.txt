cmake_minimum_required(VERSION 2.8)

find_package(MUQ_PRIVATE REQUIRED)

set(CMAKE_CXX_COMPILER ${MUQ_CXX_COMPILER})
set(CMAKE_CXX_FLAGS ${MUQ_CXX_FLAGS})

project(SimpleReadWrite)

include_directories(${MUQ_INCLUDE_DIRS})

add_executable(BlockOperations BlockOperations.cpp)
target_link_libraries(BlockOperations ${MUQ_LIBRARIES} ${MUQ_LINK_LIBRARIES})
