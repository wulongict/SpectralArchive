cmake_minimum_required(VERSION 3.2)
project(libmsms)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -std=c++11 -fPIC -fpermissive  -Werror=return-type ")
#message(MSMS the CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS})

set(SOURCE_FILES
        CArchiveSearchReply.cpp
        CThreadsPool.cpp
        ConcreteQueries.cpp
        CAnnotationDB.cpp
        CFragmentIndex.cpp
        DatabaseManager.cpp
        GraphKNN.cpp
        ICMzFile.cpp
        ICQuery.cpp
        ICGtInfoUpdate.cpp
        CMzFileReader.cpp
        mzXMLReader.cpp
        pathManager.cpp
        PeakList.cpp
        ProteomicsDataTypes.cpp
        SpectraST_base64.cpp
        SpectraST_cramp.cpp
        SWATHmzXMLReader.cpp
        SpectraST_ramp.cpp
        XMLFileParser.cpp
        Util.cpp
        CPSMAnnotation.cpp
        ICHistdata.cpp
        )
add_library(msms_lib_objects OBJECT ${SOURCE_FILES})

include_directories(../External/rapidxml-1.13 ../External/spdlog-1.x/include ../External/comet_source_2016012/MSToolkit/include )
#when using static library of sqlite3, it fail to compile.
#find_library(sqlite3  sqlite3)

add_library(msms SHARED $<TARGET_OBJECTS:msms_lib_objects> )
target_include_directories(msms INTERFACE ../External/spdlog-1.x/include ../External/comet_source_2016012/MSToolkit/include )
# target_link_libraries(msms z sqlite3 )
set_target_properties(msms PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(msms INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})



add_library(msms_static STATIC $<TARGET_OBJECTS:msms_lib_objects> )
target_include_directories(msms_static INTERFACE ../External/spdlog-1.x/include ../External/comet_source_2016012/MSToolkit/include )
target_link_libraries(msms_static  mstoolkit
 # libz.a libsqlite3.a 
 spdlog dl)
set_target_properties(msms_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(msms_static INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_dependencies(msms spdlog)
add_dependencies(msms_static spdlog)

install(TARGETS msms DESTINATION lib)
#install(FILES
#        CAnnotationDB.h        GraphKNN.h        PeakList.h
#        CArchiveSearchReply.h  ICGtInfoUpdate.h  ProteomicsDataTypes.h
#        CFragmentIndex.h       ICHistdata.h      SpectraST_base64.h
#        CMzFileReader.h        ICMzFile.h        SpectraST_ramp.h
#        ConcreteQueries.h      ICQuery.h         SWATHmzXMLReader.h
#        CThreadsPool.h         mzXMLReader.h     Util.h
#        DatabaseManager.h      pathManager.h     XMLFileParser.h
#        DESTINATION include
#        )


#---------------------------


set(SOURCE_FILES
        DatabaseManager.cpp
        CThreadsPool.cpp
        ICGtInfoUpdate.cpp
        PeakList.cpp
        ProteomicsDataTypes.cpp
        SpectraST_cramp.cpp
        SpectraST_ramp.cpp
        XMLFileParser.cpp
        Util.cpp
        CPSMAnnotation.cpp
        )

include_directories(../External/rapidxml-1.13 ../External/spdlog-1.x/include ../External/comet_source_2016012/MSToolkit/include)
#when using static library of sqlite3, it fail to compile.
#find_library(sqlite3  sqlite3)

add_library(msmsX SHARED ${SOURCE_FILES})
target_include_directories(msmsX INTERFACE ../External/spdlog-1.x/include ../External/comet_source_2016012/MSToolkit/include)
# target_link_libraries(msmsX z sqlite3 )
set_target_properties(msmsX PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(msmsX INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(msmsX_static STATIC ${SOURCE_FILES})
target_include_directories(msmsX_static INTERFACE ../External/spdlog-1.x/include ../External/comet_source_2016012/MSToolkit/include)
target_link_libraries(msmsX_static  
 # libz.a libsqlite3.a 
 spdlog dl)
set_target_properties(msmsX_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(msmsX_static INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_dependencies(msmsX spdlog)
add_dependencies(msmsX_static spdlog)

# install(TARGETS msmsX DESTINATION lib)



