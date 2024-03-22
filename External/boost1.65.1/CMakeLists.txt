
include( ExternalProject )

set(boost_zipped_file boost_1_65_1_sub.tar.gz)
# set( boost_URL "https://udomain.dl.sourceforge.net/project/boost/boost/1.65.1/boost_1_65_1.tar.gz" )
set( boost_SHA1 "b64ffb5d79da68c4bf9f8ac21b7782f7f333033d" )
set (boost_INSTALL ${CMAKE_INSTALL_PREFIX})

ExternalProject_Add( boost
        PREFIX boost
        URL ${boost_zipped_file}
        URL_HASH SHA1=${boost_SHA1}
        BUILD_IN_SOURCE 1
        # DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/External/boost1.65.1/tmp
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/External/boost1.65.1/src
        CONFIGURE_COMMAND
        ./bootstrap.sh
        --with-libraries=filesystem
        --with-libraries=system
        --with-libraries=iostreams
        --with-libraries=program_options
        --with-libraries=timer
        --with-libraries=chrono
        --prefix=<INSTALL_DIR>
       --with-python=${Python3_INCLUDE_DIRS}
        BUILD_COMMAND
#        &&
        ./b2 cxxflags=-fPIC cflags=-fPIC install
        #        link=static variant=release threading=multi runtime-link=static
        -j 30
#        echo done
        INSTALL_COMMAND ""
        # INSTALL_DIR ${boost_INSTALL}
        )


