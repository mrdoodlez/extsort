cmake_minimum_required(VERSION 2.8.4)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O2")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(extsort)

add_executable(extsort extsort.cpp cmdline.cpp)
add_executable(gendata gendata.cpp cmdline.cpp)
