
include( ExternalProject )

set( fcgi_url "http://ftp.debian.org/debian/pool/main/libf/libfcgi/libfcgi_2.4.0.orig.tar.gz" )
set( fcgi_SHA1 "940ac5e857d63756726068ac3043319289d886c9" )
set (fcgi_INSTALL ${CMAKE_INSTALL_PREFIX})

ExternalProject_Add( libfcgi
        PREFIX libfcgi
        URL ${fcgi_url}
        URL_HASH SHA1=${fcgi_SHA1}
        BUILD_IN_SOURCE 1
        DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/External/fcgi2.4.0/tmp
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/External/fcgi2.4.0/src
        CONFIGURE_COMMAND
        echo ==== <INSTALL_DIR> =====
        &&
        ./configure
        # --with-libraries=filesystem
        # --with-libraries=system
        # --with-libraries=iostreams
        # --with-libraries=program_options
        # --with-libraries=timer
        # --with-libraries=chrono
        --prefix=<INSTALL_DIR>
    #    --with-python=${Python3_INCLUDE_DIRS}
        BUILD_COMMAND
#        &&
        make 
#         ./b2 cxxflags=-fPIC cflags=-fPIC install
#         #        link=static variant=release threading=multi runtime-link=static
#         -j 30
# #        echo done
        INSTALL_COMMAND make install
        # INSTALL_DIR ${boost_INSTALL}
        )


