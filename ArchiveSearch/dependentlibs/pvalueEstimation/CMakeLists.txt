cmake_minimum_required(VERSION 3.2)
project(dependentlibs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11  -fpermissive -fopenmp -Wreturn-type  -Werror=return-type ")

# related to linear regression
#include_directories(../../../External/Eigen)



add_library(mltools SHARED EMAlgorithm.cpp Visual.cpp)
target_link_libraries(mltools ${Boost_LIBRARIES})
set_target_properties(mltools PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_dependencies(mltools boost)

add_library(mltools_static STATIC EMAlgorithm.cpp Visual.cpp)
target_link_libraries(mltools_static ${Boost_LIBRARIES})

add_dependencies(mltools_static boost)

# Everything related to tail estimation
add_library(tailestimation_static STATIC CTailEstimation.cpp  pvalueCalc.cpp )
target_link_libraries(tailestimation_static mltools_static)
add_dependencies(tailestimation_static boost)

# shared one need PIC
add_library(tailestimation SHARED CTailEstimation.cpp  pvalueCalc.cpp)
set_target_properties(tailestimation PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(tailestimation mltools)

