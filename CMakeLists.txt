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

option(ENABLE_GPU "Enable GPU based indexes, if CUDA exist" ON)
option(ENABLE_TESTS "Enable tests" ON)

if (${ENABLE_TESTS})
#	message("enable testing ...")
	# run enable_testing on root directory
	enable_testing()
endif()

#message("--- START ENABLE GPU IS " ${ENABLE_GPU} " COMPILER IS : " ${CMAKE_CUDA_COMPILER})
if(ENABLE_GPU)
	CHECK_LANGUAGE(CUDA)
	if(CMAKE_CUDA_COMPILER)
		enable_language(CUDA)
#		message("include directory when CUDA enabled: " ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
	else()
		#message("CUDA not found. Using CPU" ${CMAKE_CUDA_COMPILER})
	endif()
else()
	unset(CMAKE_CUDA_COMPILER)
	message("ENABLE GPU IS OFF " ${ENABLE_GPU} ${CMAKE_CUDA_COMPILER})
endif()

#add_subdirectory(utilities/groundtruth)
#add_subdirectory(utilities/calcrecallprecision)
#add_subdirectory(utilities/msUtility)
#add_subdirectory(utilities/readpepxml)
#add_subdirectory(utilities/extractchromatogram)
#add_subdirectory(utilities/firstcmake)
#add_subdirectory(utilities/swathnoise)
#add_subdirectory(utilities/theoreticalspectra)
#add_subdirectory(utilities/shuffleRetentionTime)
#add_subdirectory(utilities/pParse)
#add_subdirectory(utilities/LogisticRegression)

add_subdirectory(External/comet_source_2016012)
add_subdirectory(External/gtest-1.7.0)
add_subdirectory(External/liblinear-2.11)
#add_subdirectory(External/liblinear-multicore-2.42 ) # OpenMP incompatible problem. Fail to compile with old GCC versions
#add_subdirectory(External/pepnovo)
#add_subdirectory(External/ranger)
add_subdirectory(External/spdlog-1.x)
add_subdirectory(External/SpectraST)
#add_subdirectory(External/standalone-earth)
#add_subdirectory(External/svmlight)

#add_subdirectory(smallProj/NoiseFilter)
#add_subdirectory(smallProj/mars)
#add_subdirectory(smallProj/testproj)
#add_subdirectory(smallProj/SpectraSTExportor)
#add_subdirectory(smallProj/normRT)

#file(TOUCH_NOCREATE psmvalidator_main.cpp)
#add_subdirectory(psmvalidator/dependent)

add_subdirectory(librarymsms)
#add_subdirectory(unit_tests)
add_subdirectory(ArchiveSearch)

# add_subdirectory(ArchiveSearch/dependentlibs/cudaproj)


include_directories(External)
include_directories(External/spdlog-1.x/include)

set(CMAKE_CXX_FLAGS " -Wreturn-type ${CMAKE_CXX_FLAGS} -std=c++11  ")
#g++ -fsanitize=address  -fuse-ld=gold  -std=c++11 test.cpp
# this sanitize is good tool
#add_executable(mytool maintool.cpp )
#target_link_libraries(mytool -pthread log4cpp)

#add_dependencies(mzxmltomatrix msms)
#add_dependencies(fragmodel_test msbasic)
#add_dependencies(extractchrom msms)

