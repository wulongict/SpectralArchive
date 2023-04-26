cmake_minimum_required(VERSION 3.2)
project(tests)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++17 -fpermissive -fopenmp  -fPIC -Wreturn-type  -Werror=return-type ")

include_directories(../External/gtest-1.7.0/include ../librarymsms)
add_executable(msmstest msmstest.cpp)
target_link_libraries(msmstest gtest_main gtest msbasic_small  msfileramp_static)

# add_executable(new_features new_feature.cpp)

# target_compile_features(new_features PRIVATE cxx_std_17)