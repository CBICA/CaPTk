#MESSAGE( "External project - nifti" )

# if zlib is installed use it
SET( NIFTI_DEPENDENCIES )

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

#MESSAGE( STATUS "Adding nifti-${PETSC_VERSION} ...")

ExternalProject_Add(
    nifti
    SOURCE_DIR ${CMAKE_BINARY_DIR}/nifti-src
    BINARY_DIR ${CMAKE_BINARY_DIR}/nifti-build
    URL https://sourceforge.net/projects/niftilib/files/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz 
    # BUILD_COMMAND make
    INSTALL_COMMAND cmake -E echo "Skipping install step."
    CMAKE_ARGS 
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/nifti-build 
        -DCMAKE_CXX_COMPILER=mpicxx
        -DCMAKE_C_COMPILER=mpicc
        -DBUILD_SHARED_LIBS:BOOL=OFF
        -DZLIB_ROOT=${ZLIB_DIR}
        -Wno-dev
)

SET( NIFTI_DIR ${CMAKE_BINARY_DIR}/nifti-build)
#LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/OpenCV-build")
SET( ENV{CMAKE_PREFIX_PATH} "${CMAKE_PREFIX_PATH};${CMAKE_BINARY_DIR}/nifti-build" )
