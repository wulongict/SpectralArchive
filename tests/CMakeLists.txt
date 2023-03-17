cmake_minimum_required(VERSION 3.2)
project(tests)

include_directories(../External/gtest-1.7.0/include)
add_executable(msmstest msmstest.cpp)
target_link_libraries(msmstest gtest_main gtest msbasic_small  msfileramp_static)
