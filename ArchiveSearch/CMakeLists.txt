cmake_minimum_required(VERSION 3.2)
project(Spectroscape)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fpermissive -fopenmp  -fPIC -Wreturn-type  -Werror=return-type ")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)

link_directories(${CMAKE_BINARY_DIR}/Release)

# set(CMAKE_INSTALL_RPATH "${ORIGIN}/lib")

# code can be moved to other places: might rely libraries in dependentlibs
add_subdirectory(dependentlibs)
include_directories(External)

include_directories(../External/Eigen ../External/spdlog-1.x/include
		/usr/local/cuda/include
		)
set(SOURCE_FILES
		CTimerSummary.cpp
        CSpectralArchive.cpp
		FastCgiInterface.cpp
        )


if(ENABLE_GPU)
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


#The fast cgi interface
add_executable(spectroscape spectroscapeMain.cpp FileUtility.cpp dependentlibs/msbasic/commonfun.cpp dependentlibs/msbasic/CDebugMode.cpp)
target_link_libraries(spectroscape
specarchive_static archiveindex_static tailestimation_static msfileramp_static
msbasic_small_static mstoolkit spectrastlib_static
		  )
set(MY_SPECTROSCAPE_VERSION "1.1.6")
# target_compile_definitions(spectroscape PRIVATE SPECTROSCAPE_VERSION=${MY_SPECTROSCAPE_VERSION})
add_compile_definitions(SPECTROSCAPE_VERSION="${MY_SPECTROSCAPE_VERSION}")

add_dependencies(spectroscape boost ${DPGPU})

#The commandline tool for library search
#add_executable(LibraryRetrieval	LibSearchEngine.cpp	)
#target_link_libraries(LibraryRetrieval specarchive)

install(TARGETS spectroscape DESTINATION bin )
# install(TARGETS    linear  DESTINATION lib OPTIONAL)
# install(TARGETS  msfileramp msbasic_small linear  DESTINATION lib OPTIONAL)

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libopenblas-dev,libfcgi-dev,libgomp1,sqlite")
set(CPACK_PACKAGE_NAME Spectroscape_CPU)
set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local)
if(${WITH_GPU})
	set(CPACK_PACKAGE_NAME Spectroscape_GPU)
endif()

# The name of the package file to generate, not including the extension. For example, cmake-2.6.1-Linux-i686. The default value is:
# ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}
set(CPACK_PACKAGE_VERSION ${MY_SPECTROSCAPE_VERSION})
set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION})

set(CPACK_DEBIAN_PACKAGE_MAINTAINER Long WU)
set(CPACK_OUTPUT_FILE_PREFIX ${CMAKE_SOURCE_DIR}/package)
set(CPACK_THREADS 0)
SET(CPACK_GENERATOR "ZIP")  #.ZIP
SET(CPACK_GENERATOR "STGZ" "DEB") #.sh #deb
SET(CPACK_SOURCE_GENERATOR "TBZ2")

list(APPEND CPACK_SOURCE_IGNORE_FILES "/External/boost1.65.1/src")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/External/boost1.65.1/tmp" "/External/boost1.65.1/ex-*")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/.git")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/tmp")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/build")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/package")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/.vscode" ".idea" "./SpectralArchiveWeb")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/External/faiss-1.7.3/faiss_tmp*" "/External/faiss-1.7.3/html" "/External/faiss-1.7.3/xml")

list(APPEND CPACK_SOURCE_IGNORE_FILES "/cmake-build-release-cpu/.*")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/cmake-build-release-gpu/.*")
list(APPEND CPACK_SOURCE_IGNORE_FILES "/cmake-build-release/.*" "/External/SpectraST/CMakeFiles/.*" "/External/spdlog-1.x/CMakeFiles/.*")

include(CPack)






