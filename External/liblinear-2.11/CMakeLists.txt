cmake_minimum_required(VERSION 3.2)
project(liblinear)

set(CMAKE_CXX_FLAGS " -fPIC ${CMAKE_CXX_FLAGS} -std=c++11 ")

set(SOURCE_FILES
        linear.cpp
        linear.h
        tron.cpp
        tron.h
        )

set(BLAS_SOURCE_FILES
        blas/blas.h
        blas/blasp.h
        blas/dscal.c
        blas/dnrm2.c
        blas/ddot.c
        blas/daxpy.c
        )

#message(STATUS "cmake cachefile dir ${CMAKE_BINARY_DIR} -----")
#set(BLAS_STATIC_LIBRARY ${CMAKE_BINARY_DIR}/Release/libblas.a)
#message(STATUS "cmake cachefile dir ${BLAS_STATIC_LIBRARY} -----")

add_library(blas STATIC ${BLAS_SOURCE_FILES})
set_property(TARGET blas PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(train ${SOURCE_FILES} train.c )
target_link_libraries(train blas)
set_property(TARGET train PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(predict ${SOURCE_FILES} predict.c )
target_link_libraries(predict blas)
set_property(TARGET predict PROPERTY POSITION_INDEPENDENT_CODE ON)


add_library(linear SHARED ${SOURCE_FILES})
target_link_libraries(linear blas)
set_property(TARGET linear PROPERTY POSITION_INDEPENDENT_CODE ON)

add_library(linear_static STATIC ${SOURCE_FILES})
target_link_libraries(linear_static blas)
set_property(TARGET linear_static PROPERTY POSITION_INDEPENDENT_CODE ON)



# install(TARGETS train predict
#         EXPORT linear
#         DESTINATION bin)