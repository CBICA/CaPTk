#MESSAGE( "External project - petsc" )

SET( PETSC_DEPENDENCIES )

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

#MESSAGE( STATUS "Adding petsc-${PETSC_VERSION} ...")

ExternalProject_Add( 
  petsc 
  URL http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.9.4.tar.gz
  BUILD_IN_SOURCE true
  SOURCE_DIR petsc-source
  CONFIGURE_COMMAND ./configure PETSC_DIR=${CMAKE_BINARY_DIR}/petsc-source PETSC_ARCH=cxx_opt_dbl --prefix=${CMAKE_BINARY_DIR}/petsc-build/cxx_opt_dbl --with-cc=mpicc COPTFLAGS='-O3' --with-cxx=mpicxx --download-f2cblaslapack CXXOPTFLAGS='-O3' --with-ssl=0 --with-debugging=0 --with-64-bit-indices --with-shared=0 --with-x=0 --with-fc=0
  BUILD_COMMAND make PETSC_DIR=${CMAKE_BINARY_DIR}/petsc-source PETSC_ARCH=cxx_opt_dbl
  # INSTALL_COMMAND make PETSC_DIR=${CMAKE_BINARY_DIR}/petsc-source PETSC_ARCH=cxx_opt_dbl install
)

SET( PETSC_DIR ${CMAKE_BINARY_DIR}/petsc-source)
SET( PETSC_BUILD_DIR ${CMAKE_BINARY_DIR}/pets-build/cxx_opt_dbl)
#LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}/OpenCV-build")
SET( ENV{CMAKE_PREFIX_PATH} "${CMAKE_PREFIX_PATH};${CMAKE_BINARY_DIR}/petsc-source;${CMAKE_BINARY_DIR}/petsc-build/cxx_opt_dbl" )
