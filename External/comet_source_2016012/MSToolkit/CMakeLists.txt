cmake_minimum_required(VERSION 3.6)
project(MSToolkit)

set(HEADER_PATH include)


set(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -O3  -static -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DHAVE_EXPAT_CONFIG_H  ")

set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -O3 -static  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DHAVE_EXPAT_CONFIG_H  ")


include_directories(. include)


set(SOURCE_FILES

    src/expat-2.2.0/xmlparse.c
    src/expat-2.2.0/xmlrole.c
    src/expat-2.2.0/xmltok.c
    src/expat-2.2.0/xmltok_impl.c
    src/expat-2.2.0/xmltok_ns.c
    src/MSToolkit/MSObject.cpp
    src/MSToolkit/MSReader.cpp
    src/MSToolkit/mzMLWriter.cpp
    src/MSToolkit/pepXMLWriter.cpp
    src/MSToolkit/RAWReader.cpp
    src/MSToolkit/Spectrum.cpp
    src/mzParser/BasicChromatogram.cpp
    src/mzParser/BasicSpectrum.cpp
    src/mzParser/Czran.cpp
    src/mzParser/MSNumpress.cpp
    src/mzParser/mz5handler.cpp
    src/mzParser/mzMLReader.cpp
    src/mzParser/mzp_base64.cpp
    src/mzParser/mzParser.cpp
    src/mzParser/mzpMz5Config.cpp
    src/mzParser/mzpMz5Structs.cpp
    src/mzParser/PWIZface.cpp
    src/mzParser/RAMPface.cpp
    src/mzParser/saxhandler.cpp
    src/mzParser/saxmzmlhandler.cpp
    src/mzParser/saxmzxmlhandler.cpp
    src/sqlite-3.7.7.1/sqlite3.c
#    src/sqlite-3.7.7.1/sqlite3.cpp
    src/zlib-1.2.8/adler32.c
    src/zlib-1.2.8/compress.c
    src/zlib-1.2.8/crc32.c
    src/zlib-1.2.8/deflate.c
    src/zlib-1.2.8/gzclose.c
    src/zlib-1.2.8/gzlib.c
    src/zlib-1.2.8/gzread.c
    src/zlib-1.2.8/gzwrite.c
    src/zlib-1.2.8/infback.c
    src/zlib-1.2.8/inffast.c
    src/zlib-1.2.8/inflate.c
    src/zlib-1.2.8/inftrees.c
    src/zlib-1.2.8/trees.c
    src/zlib-1.2.8/uncompr.c
    src/zlib-1.2.8/zutil.c)

#add_executable(MSToolkit ${SOURCE_FILES})
add_library(mstoolkit STATIC ${SOURCE_FILES})
target_link_libraries(mstoolkit m -pthread dl)
set_target_properties(mstoolkit PROPERTIES POSITION_INDEPENDENT_CODE ON)