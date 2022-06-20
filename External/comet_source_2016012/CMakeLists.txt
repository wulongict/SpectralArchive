cmake_minimum_required(VERSION 3.6)
project(comet_source_2016012)
set(MSTOOLKIT MSToolkit)
set(COMETSEARCH CometSearch)

set(CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -fpermissive -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -I${MSTOOLKIT}/include -I${COMETSEARCH} ")
include_directories(${MSTOOLKIT}/include)
include_directories(${COMETSEARCH})
#link_directories(./MSToolkit ./CometSearch)
link_directories(${CMAKE_BINARY_DIR}/Release)

add_subdirectory(MSToolkit)
add_subdirectory(CometSearch)

set(SOURCE_FILES
    Comet.cpp
        )

add_executable(comet ${SOURCE_FILES})
target_link_libraries(comet cometsearch mstoolkit m -pthread dl)
add_dependencies(comet cometsearch mstoolkit)


