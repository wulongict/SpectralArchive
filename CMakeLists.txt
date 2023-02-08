cmake_minimum_required(VERSION 3.20)
project(mytool LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

add_subdirectory(External/boost1.65.1)
ExternalProject_Get_Property(boost install_dir )
set( boost_INCLUDE_DIR ${install_dir}/include )
set( boost_LIB_DIR ${install_dir}/lib )
set(BOOST_ROOT ${install_dir})

set(Boost_LIBRARIES
        ${boost_LIB_DIR}/libboost_filesystem.a
        ${boost_LIB_DIR}/libboost_iostreams.a
        ${boost_LIB_DIR}/libboost_program_options.a
        ${boost_LIB_DIR}/libboost_timer.a
        ${boost_LIB_DIR}/libboost_chrono.a 
        ${boost_LIB_DIR}/libboost_system.a
        )

include_directories(${boost_INCLUDE_DIR})

set(Boost_NO_SYSTEM_PATHS TRUE)
set(Boost_NO_BOOST_CMAKE TRUE)

set(Boost_FOUND TURE)

# enable GPU 
include(CheckLanguage)

set(DPGPU "")

option(ENABLE_TESTS "Enable tests" ON)

if (${ENABLE_TESTS})
	enable_testing()
endif()

if (${WITH_GPU})
        CHECK_LANGUAGE(CUDA)
endif()

if(CMAKE_CUDA_COMPILER)
	enable_language(CUDA)
	option(ENABLE_GPU "Enable GPU based indexes, if CUDA exist" ON)
	option(FAISS_ENABLE_GPU "Enable support for GPU indexes." ON)
else()
	option(FAISS_ENABLE_GPU "Enable support for GPU indexes." OFF)
	option(ENABLE_GPU "Enable GPU based indexes, if CUDA exist" OFF)
endif()


add_subdirectory(External/comet_source_2016012)
# add_subdirectory(External/gtest-1.7.0)
add_subdirectory(External/liblinear-2.11)
add_subdirectory(External/spdlog-1.x)
add_subdirectory(External/SpectraST)

option(BUILD_TESTING "build testing for faiss library" OFF)

add_subdirectory(External/faiss-1.7.1)
add_subdirectory(librarymsms)
add_subdirectory(External/gtest-1.7.0)
add_subdirectory(tests)

add_subdirectory(ArchiveSearch)

include_directories(External)
include_directories(External/spdlog-1.x/include)

set(CMAKE_CXX_FLAGS " -Wreturn-type ${CMAKE_CXX_FLAGS} -std=c++11  ")


