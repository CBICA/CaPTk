#MESSAGE( "External project - fftw" )

SET( FFTW_DEPENDENCIES )

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

#MESSAGE( STATUS "Adding fftw-${PETSC_VERSION} ...")

set(FFTW_OPTIONS --enable-sse2 MAKEINFO=missing )
if (ENABLE_OMP)
    set(FFTW_OPTIONS ${FFTW_OPTIONS} --enable-threads --enable-openmp)
endif()
if(ENABLE_AVX)
    set(FFTW_OPTIONS ${FFTW_OPTIONS} --enable-avx)
endif()
ExternalProject_Add(
    fftw-dbl
    BUILD_IN_SOURCE 1
    URL "http://www.fftw.org/fftw-3.3.6-pl2.tar.gz"
    SOURCE_DIR ${CMAKE_BINARY_DIR}/fftw-src
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/fftw-build ${FFTW_OPTIONS} CFLAGS=-O3
    BUILD_COMMAND make 
    INSTALL_COMMAND make install
)
ExternalProject_Add(
    fftw-sgl
    BUILD_IN_SOURCE 1
    URL "http://www.fftw.org/fftw-3.3.6-pl2.tar.gz"
    SOURCE_DIR ${CMAKE_BINARY_DIR}/fftw-src
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/fftw-build ${FFTW_OPTIONS} CFLAGS=-O3 --enable-float
    BUILD_COMMAND make 
    INSTALL_COMMAND make install
)

SET( FFTW_DIR ${CMAKE_BINARY_DIR}/fftw-build)
#LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/OpenCV-build")
SET( ENV{CMAKE_PREFIX_PATH} "${CMAKE_PREFIX_PATH};${CMAKE_BINARY_DIR}/fftw-build" )
