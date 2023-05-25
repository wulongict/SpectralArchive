cmake_minimum_required(VERSION 3.2)
project(msbasicproj)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wreturn-type  -fPIC -fpermissive -Werror=return-type ")

set(SOURCE_FILES
        CDebugMode.cpp
        commonfun.cpp

        )

add_library(msbasic_obj OBJECT ${SOURCE_FILES})
add_dependencies(msbasic_obj boost)

include_directories(../../../External/spdlog-1.x/include)
link_directories(${CMAKE_BINARY_DIR}/Release)

include_directories(${Boost_INCLUDE_DIRS} ../../../External/liblinear-2.11 ./)

add_library(msbasic_small SHARED $<TARGET_OBJECTS:msbasic_obj> )
set_target_properties(msbasic_small PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(msbasic_small ${Boost_LIBRARIES} pthread
        msms
        spectrastlib
        linear)
add_dependencies(msbasic_small boost)

add_library(msbasic_small_static STATIC $<TARGET_OBJECTS:msbasic_obj> )
target_link_libraries(msbasic_small_static pthread
#                msms
#        spectrastlib_static
#        linear_static
        )
#add_dependencies(msbasic_small_static boost)


