cmake_minimum_required(VERSION 3.2)
project(ArchiveSearch)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fpermissive -fopenmp  -fPIC -Wreturn-type  -Werror=return-type ")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)

#message("the path of cmake binary: " ${CMAKE_BINARY_DIR})

link_directories(${CMAKE_BINARY_DIR}/Release)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# code can be moved to other places: might rely libraries in dependentlibs
#add_subdirectory(independent_code)
add_subdirectory(dependentlibs)
include_directories(External)

include_directories(../External/Eigen ../External/spdlog-1.x/include
		/usr/local/cuda/include
		)
set(SOURCE_FILES
        CSpectralArchive.cpp
		FastCgiInterface.cpp
        )


if(CMAKE_CUDA_COMPILER)
	add_subdirectory(dependentlibs/cudaproj)
	add_compile_definitions(__CUDA__)
	set(DPGPU dpgpu)
	install(TARGETS dpgpu DESTINATION lib OPTIONAL)

endif()

add_library(specarchive SHARED ${SOURCE_FILES} )
set_target_properties(specarchive PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(specarchive fcgi fcgi++ msfileramp ${DPGPU}  archiveindex tailestimation)

add_library(specarchive_static STATIC ${SOURCE_FILES} )
target_link_libraries(specarchive_static libfcgi.a libfcgi++.a  ${DPGPU} msms   pthread)
add_dependencies(specarchive_static boost)



#The socket based method
#add_executable(similarity  spectral_similarity_gpu.cpp  CSocketServer.cpp)
#target_link_libraries(similarity specarchive)

#The fast cgi interface
add_executable(fastcgi_similarity.fcgi  similarity_gpu_fcgi.cpp )
target_link_libraries(fastcgi_similarity.fcgi    
specarchive_static archiveindex_static tailestimation_static msfileramp_static
msbasic_static faiss ${DPGPU}
    ${Boost_LIBRARIES} )
add_dependencies(fastcgi_similarity.fcgi boost ${DPGPU})

#The commandline tool for library search
#add_executable(LibraryRetrieval	LibSearchEngine.cpp	)
#target_link_libraries(LibraryRetrieval specarchive)

install(TARGETS fastcgi_similarity.fcgi DESTINATION bin OPTIONAL)
install(TARGETS  msfileramp msbasic linear  DESTINATION lib OPTIONAL)
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 0)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_NAME archiveSearch)
set(CPACK_THREADS 0)
SET(CPACK_GENERATOR "ZIP")  #.ZIP
SET(CPACK_GENERATOR "STGZ") #.sh
SET(CPACK_SOURCE_GENERATOR "TBZ2")

list(APPEND CPACK_SOURCE_IGNORE_FILES "/External/boost1.65.1/src")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/External/boost1.65.1/tmp")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/.git")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/tmp")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/build")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/cmake-build-release/.*")

include(CPack)






