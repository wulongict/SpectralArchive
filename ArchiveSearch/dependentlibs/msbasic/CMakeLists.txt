cmake_minimum_required(VERSION 3.2)
project(msbasicproj)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wreturn-type  -fPIC -fpermissive -Werror=return-type ")

set(SOURCE_FILES
        BasePSMFeature.cpp
        CDebugMode.cpp
        CFragScore.cpp
        classifier.cpp
        commonfun.cpp
        ConcretePSMFeatures.cpp
        CPeakPair.cpp
        CPeakPairsImporter.cpp
        FragmentationScore.cpp
        gnuplot_functions.cpp
        liblinearWrapper.cpp
        MGFReader.cpp
        PeakInfo.cpp
        PeakPairScore.cpp
        peak_annotation.cpp
        PSMScoreCalculator.cpp
        CFlow.cpp
        )

add_library(msbasic_obj OBJECT ${SOURCE_FILES})
add_dependencies(msbasic_obj boost)

include_directories(../../../External/spdlog-1.x/include)
link_directories(${CMAKE_BINARY_DIR}/Release)

include_directories(${Boost_INCLUDE_DIRS} ../../../External/liblinear-2.11 ./)

add_library(msbasic SHARED $<TARGET_OBJECTS:msbasic_obj> )
set_target_properties(msbasic PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(msbasic ${Boost_LIBRARIES} pthread
        msms
        spectrastlib
        linear)
add_dependencies(msbasic boost)

add_library(msbasic_static STATIC $<TARGET_OBJECTS:msbasic_obj> )
target_link_libraries(msbasic_static pthread
                msms
        spectrastlib_static
        linear_static
        )
add_dependencies(msbasic_static boost)


