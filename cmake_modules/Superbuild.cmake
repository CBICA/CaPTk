## The Superbuild script is used to automatically download and build project dependencies. 

# Using GIT to download third party libraries. An SVN/BitBucket URL will also work the same way
FIND_PACKAGE( Git REQUIRED )

OPTION( USE_GIT_PROTOCOL "If behind a firewall turn this off to use https instead." OFF )

SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/qt/5.11.2" )
SET( FILENAME_TO_EXTRACT "qt" )
SET( FILE_TO_EXTRACT "${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip" )
SET( QT_EXTRACTED_DIR "${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}" )

IF( NOT EXISTS "${QT_EXTRACTED_DIR}" )
  FILE(MAKE_DIRECTORY "${QT_EXTRACTED_DIR}" )
ENDIF()

IF(WIN32)
  SET( DOWNLOAD_LINK "${DOWNLOAD_LINK}/windows.zip" )
ELSEIF(APPLE)
  SET( DOWNLOAD_LINK "${DOWNLOAD_LINK}/macos.zip" )
ELSE()
  SET( DOWNLOAD_LINK "${DOWNLOAD_LINK}/linux.zip" )
ENDIF()

IF( NOT EXISTS "${FILE_TO_EXTRACT}" )
  MESSAGE( STATUS "Downloading pre-compiled Qt with open source license (see Qt site for more details)" )
  FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000 STATUS STATUS_CODE SHOW_PROGRESS)
  IF(NOT STATUS_CODE EQUAL 0)
    MESSAGE(FATAL_ERROR "Failed to download Precompiled packages. Status=${STATUS_CODE}")
  ENDIF()

ENDIF()

IF( EXISTS "${FILE_TO_EXTRACT}" )

  IF( NOT EXISTS "${QT_EXTRACTED_DIR}/5.11.2/lib/cmake/Qt5" )
    MESSAGE( STATUS "Extracting pre-compiled Qt binaries" )
    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
      WORKING_DIRECTORY ${QT_EXTRACTED_DIR}
      RESULT_VARIABLE RESULT_CODE
    )

    IF(NOT RESULT_CODE EQUAL 0)
      MESSAGE( WARNING "Extracting the pre-compiled Qt binaries failed" )
    ENDIF()
  ENDIF()

  LIST(APPEND CMAKE_PREFIX_PATH "${QT_EXTRACTED_DIR}/5.11.2/lib/cmake/Qt5")
  LIST(APPEND CMAKE_PROGRAM_PATH  "${QT_EXTRACTED_DIR}/5.11.2/bin")
  
  FIND_PACKAGE(Qt5 COMPONENTS Core Gui Svg Widgets WebView WebEngine WebEngineCore)

ENDIF()
  
LINK_DIRECTORIES(${QT_LIBRARY_DIR})
SET( Qt5_DIR ${Qt5_DIR} CACHE STRING "Qt5_DIR Path for use in other builds" FORCE )
LIST(APPEND CMAKE_PREFIX_PATH "${Qt5_DIR}")

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

#MESSAGE( STATUS "Adding DCMTK-3.6.3 ...")
#INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-DCMTK.cmake )

#MESSAGE( STATUS "Adding YAML-CPP 0.6.2 ...")
#INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-yaml-cpp.cmake )

MESSAGE( STATUS "Adding EIGEN-3.3.4 ...")
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-Eigen.cmake )

MESSAGE( STATUS "Adding OpenCV-3.4.1 ...")
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-OpenCV.cmake )

MESSAGE( STATUS "Adding VTK-8.1.0 ...")
INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-VTK.cmake )

#MESSAGE( STATUS "Adding ITK-4.13.0 ...")
#INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/External-ITK.cmake )
