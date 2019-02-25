## The Superbuild script is used to automatically download and build project dependencies. 

# Using GIT to download third party libraries. An SVN/BitBucket URL will also work the same way
FIND_PACKAGE( Git REQUIRED )

OPTION( USE_GIT_PROTOCOL "If behind a firewall turn this off to use https instead." OFF )

#LIST(APPEND CMAKE_PROGRAM_PATH  "c:/MyProject/Tools/mingw/bin/" ...)

SET(git_protocol "git")
IF(NOT USE_GIT_PROTOCOL)
	SET(git_protocol "https")
ENDIF()

SET( 
  CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake_modules
  ${CMAKE_MODULE_PATH}
)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

IF(MSVC)
  SET( CMAKE_BUILD_TYPE "Debug;Release")
ELSEIF(UNIX)
  SET( CMAKE_BUILD_TYPE "Release")
ENDIF()

INCLUDE( ExternalProject )

## Compute -G arg for configuring external projects with the same CMake generator:
#IF(CMAKE_EXTRA_GENERATOR)
#	SET(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
#ELSE()
#	SET(gen "${CMAKE_GENERATOR}" )
#ENDIF()

#INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-DCMTK.cmake )

INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-Eigen.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-VTK.cmake )
#INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-OpenCV_Contrib.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-OpenCV.cmake )
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-ITK.cmake )