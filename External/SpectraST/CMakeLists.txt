cmake_minimum_required(VERSION 3.2)
project(Myspectrast)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fpermissive ")
set(CMAKE_CXX_LINK_FLAGS "-lz")

## we need gsl of 1.16 version not higher!!
## will fix later, but now let's check the version of GSL first
#
#include(FindGSL)
#message(STATUS "GSL VERSION ${GSL_VERSION}")
#set(GSL_VER_REQ 1.16)
#if(GSL_FOUND )
#    message(STATUS "GSL VERSION ${GSL_VERSION}")
#    if(${GSL_VERSION} STREQUAL ${GSL_VER_REQ})
#        set(COMPILE_SPECTRAST YES)
##        message(STATUS SpectraST: found GSL ${GSL_INCLUDE_DIRS} ${GSL_LIBRARIES} )
#    else()
#        message(FATAL_ERROR "SpectraST: GSL version should be ${GSL_VER_REQ}, while we found version ${GSL_VERSION}")
#    endif()
#else()
#    message(STATUS "SpectraST: GSL NOT found! Will not compile SpectraST")
#endif()

set(RAMPSOURCEFILES
        SpectraST_base64.cpp
        SpectraST_cramp.cpp
        SpectraST_ramp.cpp
        SpectraST_util.cpp
        )

#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE Release)
#endif()

set(HPPFILES
        FileUtils.hpp
        Peptide.hpp
        Predicate.hpp
        ProgressCount.hpp
        SpectraSTCandidate.hpp
        SpectraSTConstants.hpp
        SpectraST_cramp.hpp
        SpectraSTCreateParams.hpp
        SpectraSTDenoiser.hpp
        SpectraSTDtaSearchTask.hpp
        SpectraSTFastaFileHandler.hpp
        SpectraSTFileList.hpp
        SpectraSTHtmlSearchOutput.hpp
        SpectraSTLibEntry.hpp
        SpectraSTLib.hpp
        SpectraSTLibImporter.hpp
        SpectraSTLibIndex.hpp
        SpectraSTLog.hpp
        SpectraSTMgfSearchTask.hpp
        SpectraSTMs2LibImporter.hpp
        SpectraSTMspLibImporter.hpp
        SpectraSTMspSearchTask.hpp
        SpectraSTMzLibIndex.hpp
        SpectraSTMzXMLLibImporter.hpp
        SpectraSTMzXMLSearchTask.hpp
        SpectraSTPeakList.hpp
        SpectraSTPeptideLibIndex.hpp
        SpectraSTPepXMLLibImporter.hpp
        SpectraSTPepXMLSearchOutput.hpp
        SpectraSTQuery.hpp
        SpectraSTReplicates.hpp
        SpectraSTSearch.hpp
        SpectraSTSearchOutput.hpp
        SpectraSTSearchParams.hpp
        SpectraSTSearchTask.hpp
        SpectraSTSearchTaskStats.hpp
        SpectraSTSimScores.hpp
        SpectraSTSpLibImporter.hpp
        SpectraSTTsvLibImporter.hpp
        SpectraSTTxtSearchOutput.hpp
        SpectraSTXHunterLibImporter.hpp
        SpectraSTXlsSearchOutput.hpp
        XMLWalker.hpp
        )

set(SOURCEFILES
        #cppnames
        FileUtils.cpp
        FileUtils.hpp
        #Lib2HTML.cpp
        Peptide.cpp
        Peptide.hpp
        #plotspectrast.cpp
        Predicate.cpp
        Predicate.hpp
        ProgressCount.cpp
        ProgressCount.hpp
        SpectraST_base64.cpp
        SpectraSTCandidate.cpp
        SpectraSTCandidate.hpp
        SpectraSTConstants.hpp
        SpectraST_cramp.cpp
        SpectraST_cramp.hpp
        SpectraSTCreateParams.cpp
        SpectraSTCreateParams.hpp
        SpectraSTDenoiser.cpp
        SpectraSTDenoiser.hpp
        SpectraSTDtaSearchTask.cpp
        SpectraSTDtaSearchTask.hpp
        SpectraSTFastaFileHandler.cpp
        SpectraSTFastaFileHandler.hpp
        SpectraSTFileList.cpp
        SpectraSTFileList.hpp
        SpectraSTHtmlSearchOutput.cpp
        SpectraSTHtmlSearchOutput.hpp
        SpectraST_kwset.cpp
        SpectraSTLib.cpp
        SpectraSTLibEntry.cpp
        SpectraSTLibEntry.hpp
        SpectraSTLib.hpp
        SpectraSTLibImporter.cpp
        SpectraSTLibImporter.hpp
        SpectraSTLibIndex.cpp
        SpectraSTLibIndex.hpp
        SpectraSTLog.cpp
        SpectraSTLog.hpp
        #SpectraSTMain.cpp
        SpectraSTMgfSearchTask.cpp
        SpectraSTMgfSearchTask.hpp
        SpectraSTMs2LibImporter.cpp
        SpectraSTMs2LibImporter.hpp
        SpectraSTMspLibImporter.cpp
        SpectraSTMspLibImporter.hpp
        SpectraSTMspSearchTask.cpp
        SpectraSTMspSearchTask.hpp
        SpectraSTMzLibIndex.cpp
        SpectraSTMzLibIndex.hpp
        SpectraSTMzXMLLibImporter.cpp
        SpectraSTMzXMLLibImporter.hpp
        SpectraSTMzXMLSearchTask.cpp
        SpectraSTMzXMLSearchTask.hpp
        SpectraST_obstack.cpp
        SpectraSTPeakList.cpp
        SpectraSTPeakList.hpp
        SpectraSTPeptideLibIndex.cpp
        SpectraSTPeptideLibIndex.hpp
        SpectraSTPepXMLLibImporter.cpp
        SpectraSTPepXMLLibImporter.hpp
        SpectraSTPepXMLSearchOutput.cpp
        SpectraSTPepXMLSearchOutput.hpp
        SpectraSTQuery.cpp
        SpectraSTQuery.hpp
        SpectraST_ramp.cpp
        SpectraSTReplicates.cpp
        SpectraSTReplicates.hpp
        SpectraSTSearch.cpp
        SpectraSTSearch.hpp
        SpectraSTSearchOutput.cpp
        SpectraSTSearchOutput.hpp
        SpectraSTSearchParams.cpp
        SpectraSTSearchParams.hpp
        SpectraSTSearchTask.cpp
        SpectraSTSearchTask.hpp
        SpectraSTSearchTaskStats.cpp
        SpectraSTSearchTaskStats.hpp
        SpectraSTSimScores.cpp
        SpectraSTSimScores.hpp
        SpectraSTSpLibImporter.cpp
        SpectraSTSpLibImporter.hpp
        SpectraSTTsvLibImporter.cpp
        SpectraSTTsvLibImporter.hpp
        SpectraSTTxtSearchOutput.cpp
        SpectraSTTxtSearchOutput.hpp
        SpectraST_util.cpp
        SpectraSTXHunterLibImporter.cpp
        SpectraSTXHunterLibImporter.hpp
        SpectraSTXlsSearchOutput.cpp
        SpectraSTXlsSearchOutput.hpp
        XMLWalker.cpp
        XMLWalker.hpp
        )

set(PLOTSPECTRASTFILES
        plotspectrast.cpp
        SpectraSTLibEntry.cpp
        SpectraSTQuery.cpp
        SpectraSTPeakList.cpp
        SpectraSTDenoiser.cpp
        SpectraSTLog.cpp
        FileUtils.cpp
        Peptide.cpp
        ${RAMPSOURCEFILES}
        )

set(LIBHTMLSOURCEFILES
        SpectraSTPeptideLibIndex.cpp
        SpectraSTMzLibIndex.cpp
        SpectraSTLibIndex.cpp
        SpectraSTLibEntry.cpp
        SpectraSTPeakList.cpp
        SpectraSTDenoiser.cpp
        SpectraSTLog.cpp
        FileUtils.cpp
        Peptide.cpp
        ${RAMPSOURCEFILES}
        ${HPPFILES}
        )

#set( PROJECT_LINK_LIBS libmsms.so )

if(COMPILE_SPECTRAST)
#    MESSAGE(STATUS "targets: spectrast, plogspectrast Lib2HTML and spectrastlib")
    add_executable(spectrast ${SOURCEFILES} SpectraSTMain.cpp)
    target_link_libraries(spectrast  z pthread)
    target_compile_options(spectrast PRIVATE -D__LGPL__ -DSTANDALONE_LINUX)

    add_executable(plotspectrast ${PLOTSPECTRASTFILES} plotspectrast.cpp)
    target_link_libraries(plotspectrast  z pthread)
    target_compile_options(plotspectrast PRIVATE -D__LGPL__ -DSTANDALONE_LINUX)

    add_executable(plotspectrast.cgi ${PLOTSPECTRASTFILES} plotspectrast.cpp)
    target_link_libraries(plotspectrast.cgi z pthread)
    target_compile_options(plotspectrast.cgi PRIVATE -D__LGPL__  -DRUN_AS_CGI  -DSTANDALONE_LINUX)


    add_executable(Lib2HTML ${LIBHTMLSOURCEFILES} Lib2HTML.cpp)
    target_link_libraries(Lib2HTML  z pthread)
    target_compile_options(Lib2HTML  PRIVATE -D__LGPL__  -DSTANDALONE_LINUX)


    add_library(spectrastlib SHARED ${SOURCEFILES} SpectraSTMain.cpp)
#    target_link_libraries(spectrastlib gsl gslcblas z pthread)
    target_link_libraries(spectrastlib  z pthread) # do not use GSL
    set_target_properties(spectrastlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_compile_options(spectrastlib PRIVATE -D__LGPL__ -DLIB_BUILD -DSTANDALONE_LINUX)


#    target_link_options(spectrastlib PRIVATE -static)


    # install(TARGETS spectrastlib DESTINATION lib)

    add_library(spectrastlib_static STATIC ${SOURCEFILES} SpectraSTMain.cpp)
#    target_link_libraries(spectrastlib_static libgsl.a libgslcblas.a libz.a pthread dl)
    target_link_libraries(spectrastlib_static libz.a pthread dl)  # do not use GSL
#    target_link_options(spectrastlib_static PRIVATE -static)
    target_compile_options(spectrastlib_static PUBLIC -D__LGPL__ -DLIB_BUILD -DSTANDALONE_LINUX)

endif ()





