## The Superbuild script is used to automatically download and build project dependencies. 

# Using GIT to download third party libraries. An SVN/BitBucket URL will also work the same way
FIND_PACKAGE( Git REQUIRED )

OPTION( USE_GIT_PROTOCOL "If behind a firewall turn this off to use https instead." OFF )

IF(MSVC)
  SET( CMAKE_CONFIGURATION_TYPES "Debug;Release")
ELSEIF(UNIX)
  SET( CMAKE_CONFIGURATION_TYPES "Release")
  SET( CMAKE_BUILD_TYPE "Release" )
ENDIF()

INCLUDE( ExternalProject )

# check build path lenght for windows and give a warning if greater than 15
IF( WIN32 )
  STRING( LENGTH ${PROJECT_BINARY_DIR} BUILD_PATH_LENGTH )
  IF( ${BUILD_PATH_LENGTH} GREATER 15 )
    MESSAGE( WARNING "WARNING: The Superbuild path is greater than 15; it is HIGHLY recommended to make this shorter so that ITK library linkage will succeed" )
  ENDIF()
ENDIF()

## Compute -G arg for configuring external projects with the same CMake generator:
#IF(CMAKE_EXTRA_GENERATOR)
#	SET(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}") 
#ELSE()
#	SET(gen "${CMAKE_GENERATOR}" )
#ENDIF()

INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-DCMTK.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-Qt.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-VTK.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-Eigen.cmake )
#INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-OpenCV_Contrib.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-OpenCV.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-ITK.cmake )
IF( BUILD_CLAIRE )
  INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-PETSC.cmake )
  INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-fftw.cmake )
  INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-accfft.cmake )
  INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-zlib.cmake )
  INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-nifti.cmake )
ENDIF()
# INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-SEM.cmake )
