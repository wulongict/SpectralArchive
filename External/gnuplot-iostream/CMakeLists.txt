cmake_minimum_required(VERSION 3.2)
project(example_misc)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fpermissive -fopenmp ")
#set(CMAKE_CXX_LINK_FLAGS "-lz")


#link_directories(${CMAKE_BINARY_DIR}/Release)


set(SOURCE_FILES example-misc.cc )

add_executable(example_misc ${SOURCE_FILES})
target_link_libraries(example_misc boost_iostreams boost_system boost_filesystem)
target_link_directories(example_misc ${CMAKE_BINARY_DIR}/Release)
#http://stackoverflow.com/questions/10075810/gnuplot-iostream-not-linking-to-boost

