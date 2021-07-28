CMAKE_MINIMUM_REQUIRED(VERSION 3.7.2)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

SET( PROJECT_NAME Utilities )

CAPTK_ADD_PROJECT( ${PROJECT_NAME} 1.6.2 )

ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/ cbica_toolkit ) # takes care of find_itk

CAPTK_INITIAL_SETUP()

## dcmqi stuff starts

FIND_PACKAGE( DCMTK REQUIRED )

# update submodule
IF( NOT EXISTS "${PROJECT_SOURCE_DIR}/dcmqi/CMakeLists.txt" )
  EXECUTE_PROCESS(COMMAND git submodule update --init ${PROJECT_SOURCE_DIR}/dcmqi
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RESULT_VARIABLE git_result
  OUTPUT_VARIABLE git_ver)
ENDIF()

SET( CACHED_INCLUDE_DIRS
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/dcmqi/include/
  ${PROJECT_SOURCE_DIR}/dcmqi/include/dcmqi
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/json
  ${DCMTK_INCLUDE_DIRS}
)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/dcmqi/include/dcmqi)

set(HDRS
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/json/json.h
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/json/json-forwards.h
  ${INCLUDE_DIR}/JSONMetaInformationHandlerBase.h
  ${INCLUDE_DIR}/JSONParametricMapMetaInformationHandler.h
  ${INCLUDE_DIR}/JSONSegmentationMetaInformationHandler.h
  ${INCLUDE_DIR}/preproc.h
  ${INCLUDE_DIR}/QIICRConstants.h
  ${INCLUDE_DIR}/QIICRUIDs.h
  ${INCLUDE_DIR}/ConverterBase.h
  ${INCLUDE_DIR}/Exceptions.h
  ${INCLUDE_DIR}/framesorter.h
  ${INCLUDE_DIR}/ImageSEGConverter.h
  ${INCLUDE_DIR}/ParaMapConverter
  ${INCLUDE_DIR}/Helper.h
  ${INCLUDE_DIR}/SegmentAttributes.h
  ${INCLUDE_DIR}/TID1500Reader.h
  ${PROJECT_SOURCE_DIR}/dcmqi/include/cbicaDCMQIWrapper.h
  )

set(SRCS
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/JSONMetaInformationHandlerBase.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/JSONParametricMapMetaInformationHandler.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/JSONSegmentationMetaInformationHandler.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/jsoncpp.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/ConverterBase.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/ImageSEGConverter.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/ParaMapConverter.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/Helper.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/SegmentAttributes.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/TID1500Reader.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/cbicaDCMQIWrapper.cpp
  )
  
INCLUDE_DIRECTORIES( ${CACHED_INCLUDE_DIRS} )

set(lib_name dcmqi)

add_library(${lib_name} STATIC
  ${HDRS}
  ${SRCS}
  ${ADDITIONAL_SRCS}
  )
  
target_link_libraries(${lib_name}
  ${DCMTK_LIBRARIES}
  ${ITK_LIBRARIES}
)
## dcmqi stuff ends

CAPTK_ADD_EXECUTABLE( ${PROJECT_NAME} "${APPLICATION_SOURCES};${APPLICATION_DEPENDS}" "${LIBNAME_CBICATK};${lib_name}" )

ADD_DEPENDENCIES( ${PROJECT_NAME} ${LIBNAME_CBICATK} ${lib_name} ) 

# for hausdorff
SET( PLATFORM_STRING "" )

IF(WIN32)
  SET( PLATFORM_STRING "windows" )
ELSEIF(APPLE)
  SET( PLATFORM_STRING "macos" )
ELSE()
  SET( PLATFORM_STRING "linux" )
ENDIF()

SET( DOWNLOAD_LINK "https://captk.projects.nitrc.org/Hausdorff95_${PLATFORM_STRING}.zip" )
SET( FILENAME_TO_EXTRACT "hausdorff95_${PLATFORM_STRING}")
SET( FILE_TO_EXTRACT ${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip)

SET( DOWNLOADED_APPS_DIR ${PROJECT_BINARY_DIR}/hausdorff95 )

IF( NOT EXISTS "${DOWNLOADED_APPS_DIR}" )

  IF( NOT EXISTS "${FILE_TO_EXTRACT}" )

    # do not re-download if the LFS fetch worked
    IF(NOT EXISTS ${FILE_TO_EXTRACT})
      # download exe from url
      MESSAGE( STATUS "Downloading pre-compiled hausdorff" )
      FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
      IF(NOT STATUS_CODE EQUAL 0)
        SET( DOWNLOAD_LINK "https://github.com/CBICA/precompiledBinaries/raw/master/${PLATFORM_STRING}/Hausdorff95_${PLATFORM_STRING}.zip" )
        FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
        IF(NOT STATUS_CODE EQUAL 0)
          MESSAGE(FATAL_ERROR "Failed to download Precompiled hausdorff. Status=${STATUS_CODE}")
        ENDIF()
      ENDIF()
    ENDIF()
  ENDIF()

  FILE(MAKE_DIRECTORY ${DOWNLOADED_APPS_DIR})

  MESSAGE( STATUS "Extracting pre-compiled hausdorff" )
  IF( EXISTS "${FILE_TO_EXTRACT}" )

    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
      WORKING_DIRECTORY ${DOWNLOADED_APPS_DIR}
      RESULT_VARIABLE RESULT_CODE
    )

    IF(NOT RESULT_CODE EQUAL 0)
      MESSAGE( WARNING "Extracting the pre-compiled hausdorff failed" )
    ENDIF()

  ENDIF()
ENDIF()

FILE(GLOB HausdoffEXE "${DOWNLOADED_APPS_DIR}/*")

IF(APPLE)
  INSTALL(FILES ${HausdoffEXE}
    DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
ELSE()
  INSTALL(FILES ${HausdoffEXE}
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
ENDIF()
