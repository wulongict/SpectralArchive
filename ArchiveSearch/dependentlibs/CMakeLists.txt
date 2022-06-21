cmake_minimum_required(VERSION 3.2)
project(dependentlibs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11  -fpermissive -fopenmp -Wreturn-type  -Werror=return-type ")


set(MSTOOLKIT ../../External/comet_source_2016012/MSToolkit)


add_library(msfileramp SHARED randomSpecParser.cpp MzRTMap.cpp Clustering.cpp)
target_include_directories(msfileramp PRIVATE ${MSTOOLKIT}/include )
target_link_libraries(msfileramp mstoolkit msbasic)
target_compile_definitions(msfileramp PUBLIC STANDALONE_LINUX)
target_compile_options(msfileramp PRIVATE -DLIB_BUILD -O2 -Wall  -Wextra  -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__)

set_target_properties(msfileramp PROPERTIES POSITION_INDEPENDENT_CODE ON)


add_library(msfileramp_static STATIC randomSpecParser.cpp MzRTMap.cpp Clustering.cpp)
target_include_directories(msfileramp_static PRIVATE ${MSTOOLKIT}/include )
target_link_libraries(msfileramp_static mstoolkit   pthread)
target_compile_definitions(msfileramp_static PUBLIC STANDALONE_LINUX)
target_compile_options(msfileramp_static PRIVATE -DLIB_BUILD -O2 -Wall  -Wextra  -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__)



#add_subdirectory(../../External/comet_source_2016012)
add_subdirectory(spectralIndex)
add_subdirectory(pvalueEstimation)
add_subdirectory(msbasic)

