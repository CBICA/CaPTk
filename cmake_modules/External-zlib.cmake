#MESSAGE( "External project - zlib" )

SET( ZLIB_DEPENDENCIES )

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

#MESSAGE( STATUS "Adding zlib-${PETSC_VERSION} ...")

ExternalProject_Add(
    zlib
    BUILD_IN_SOURCE 1
    SOURCE_DIR ${CMAKE_BINARY_DIR}/zlib-src
    URL https://zlib.net/zlib-1.2.11.tar.gz 
    # BUILD_COMMAND make
    INSTALL_COMMAND cmake -E echo "Skipping install step."
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/zlib-build --static
)

SET( ZLIB_DIR ${CMAKE_BINARY_DIR}/zlib-build)
#LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/OpenCV-build")
SET( ENV{CMAKE_PREFIX_PATH} "${CMAKE_PREFIX_PATH};${CMAKE_BINARY_DIR}/zlib-build" )
