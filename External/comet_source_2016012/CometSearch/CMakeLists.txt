cmake_minimum_required(VERSION 3.6)
project(comet_source_2016012)
set(MSTOOLKIT ./MSToolkit)
set(COMETSEARCH ./CometSearch)

set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-write-strings -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -I. -I${MSTOOLKIT}/include ")
include_directories(${MSTOOLKIT}/include)
include_directories(${COMETSEARCH})




set(SOURCE_FILES
        CometInterfaces.cpp
        CometMassSpecUtils.cpp
        CometPostAnalysis.cpp
        CometPreprocess.cpp
        CometSearch.cpp
        CometSearchManager.cpp
       CometWriteOut.cpp
        CometWritePepXML.cpp
        CometWritePercolator.cpp
        CometWriteSqt.cpp
        CometWriteTxt.cpp
        Threading.cpp
        )

add_library(cometsearch STATIC ${SOURCE_FILES})
target_link_libraries(cometsearch m -pthread dl)

#target_link_libraries(comet_source_2016012 -lm -pthread -ldl)