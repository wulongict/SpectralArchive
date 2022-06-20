cmake_minimum_required(VERSION 3.2)
project(dependentlibs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11  -fpermissive -fopenmp -Wreturn-type  -Werror=return-type ")

#related to Eigen
include_directories(../../../External/Eigen ../../../External/spdlog-1.x/include)
add_library(kmeans_eigen_static STATIC CKMeans.cpp)
set_target_properties(kmeans_eigen_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(SOURCE_FILES 
        CFaissIndex.cpp
        ICIndexWrapper.cpp
        CMyIndex.cpp
        CKMeans.cpp
        CMultiIndices.cpp
        GpuResourceObj.cpp
)

add_library(archiveindex_obj OBJECT ${SOURCE_FILES})
add_dependencies(archiveindex_obj boost)
target_include_directories(archiveindex_obj PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
        )


# Everything related to Index
if(CMAKE_CUDA_COMPILER)
    add_compile_definitions(__CUDA__)
else()
    message(STATUS "Using CPU: GPU_FOUND=${GPU_FOUND}" )
endif()

add_library(archiveindex SHARED $<TARGET_OBJECTS:archiveindex_obj> )
set_target_properties(archiveindex PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(archiveindex  PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
        )
target_link_libraries(archiveindex faiss lapack)
add_dependencies(archiveindex boost)

add_library(archiveindex_static STATIC $<TARGET_OBJECTS:archiveindex_obj> )
target_link_libraries(archiveindex_static faiss lapack)
add_dependencies(archiveindex_static boost)

