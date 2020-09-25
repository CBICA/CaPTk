CMAKE_MINIMUM_REQUIRED(VERSION 3.7.2)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

SET( PROJECT_NAME FeatureExtraction )

CAPTK_ADD_PROJECT( ${PROJECT_NAME} 1.6.2 )

ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/ cbica_toolkit ) # takes care of find_itk

CAPTK_INITIAL_SETUP()

INCLUDE_DIRECTORIES( ${CACHED_INCLUDE_DIRS} )

# For new feature, add the feature definition name here; it needs to correspond to the file in "${DATA_DIR}/features"
SET( FeatureDefinitions "" )
SET( FT_META_FILES_TO_INCLUDE "")
SET( BATCHMODE "${PROJECT_SOURCE_DIR}/data/batchMode/batch_featureExtraction.csv" )
FILE( GLOB allFeatureDefinitionFiles "${PROJECT_SOURCE_DIR}/data/*.csv" )
FOREACH(featureDefinitionFile ${allFeatureDefinitionFiles})

  GET_FILENAME_COMPONENT( CURRENT_DEF_FILE ${featureDefinitionFile} NAME_WE )

  IF( APPLE )
    SET_SOURCE_FILES_PROPERTIES(${featureDefinitionFile} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/features)
  ELSE()
    INSTALL(
      FILES "${featureDefinitionFile}"
      DESTINATION data/features
    )
  ENDIF()

  SET( FeatureDefinitions "${FeatureDefinitions} ${CURRENT_DEF_FILE}" )
ENDFOREACH()

IF (APPLE)  
  SET_SOURCE_FILES_PROPERTIES(${allFeatureDefinitionFiles} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/features)

  SET_SOURCE_FILES_PROPERTIES(${BATCHMODE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/share/featureExtractionBatch)

  SET(FT_META_FILES_TO_INCLUDE ${FT_META_FILES_TO_INCLUDE} ${allFeatureDefinitionFiles} ${BATCHMODE})

ELSE()
  INSTALL( FILES ${BATCHMODE} DESTINATION share/featureExtractionBatch )
ENDIF()

SET( FeatureDefinitions ${FeatureDefinitions} CACHE STRING "Available Feature Definitions" FORCE)
ADD_DEFINITIONS(-DFeatureDefinitions="${FeatureDefinitions}")

SET( FEATURE_EXTRACTION_META_FILES_TO_INCLUDE "${FT_META_FILES_TO_INCLUDE}" CACHE STRING "macOS features extraction" FORCE)

CAPTK_ADD_EXECUTABLE( ${PROJECT_NAME} "${APPLICATION_SOURCES};${APPLICATION_DEPENDS}" ${LIBNAME_CBICATK} )

ADD_DEPENDENCIES( ${PROJECT_NAME} ${LIBNAME_CBICATK} )

OPTION( BUILD_REFERENCES "Build reference executables from CBIG pipeline" OFF )

IF( BUILD_REFERENCES )

  SET_PROPERTY( GLOBAL PROPERTY USE_FOLDERS ON )
  SET( StandAloneCLIAppsFolder "References")
  # find stand-alone CMD applications
  FILE( GLOB STANDALONE_APPS
    "${PROJECT_SOURCE_DIR}/src/reference/*.cxx"
  )

  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/src/reference/ )
    
  FOREACH(application ${STANDALONE_APPS})

    GET_FILENAME_COMPONENT( STANDALONE_APP_NAME ${application} NAME_WE )

    ADD_EXECUTABLE( 
      ${STANDALONE_APP_NAME}
      ${PROJECT_SOURCE_DIR}/src/reference/${STANDALONE_APP_NAME}.cxx
    )
    
    TARGET_LINK_LIBRARIES( 
      ${STANDALONE_APP_NAME}
      ${ITK_LIBRARIES}
      ${OpenCV_LIBRARIES}
    )
    
    SET_TARGET_PROPERTIES( ${STANDALONE_APP_NAME} PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )
    
  ENDFOREACH()
  
ENDIF()

# for collage
SET( PLATFORM_STRING "" )

IF(WIN32)
  SET( PLATFORM_STRING "windows" )
ELSEIF(APPLE)
  SET( PLATFORM_STRING "macos" )
ELSE()
  SET( PLATFORM_STRING "linux" )
ENDIF()

SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/collage_${PLATFORM_STRING}.zip" )
SET( FILENAME_TO_EXTRACT "collage_${PLATFORM_STRING}")
SET( FILE_TO_EXTRACT ${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip)

SET( DOWNLOADED_APPS_DIR ${PROJECT_BINARY_DIR}/collage )

IF( NOT EXISTS "${DOWNLOADED_APPS_DIR}" )

  IF( NOT EXISTS "${FILE_TO_EXTRACT}" )

    # do not re-download if the LFS fetch worked
    IF(NOT EXISTS ${FILE_TO_EXTRACT})
      # download exe from url
      MESSAGE( STATUS "Downloading pre-compiled collage" )
      FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
      IF(NOT STATUS_CODE EQUAL 0)
        MESSAGE(FATAL_ERROR "Failed to download Precompiled collage. Status=${STATUS_CODE}")
      ENDIF()
    ENDIF()
  ENDIF()

  FILE(MAKE_DIRECTORY ${DOWNLOADED_APPS_DIR})

  MESSAGE( STATUS "Extracting pre-compiled collage" )
  IF( EXISTS "${FILE_TO_EXTRACT}" )

    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
      WORKING_DIRECTORY ${DOWNLOADED_APPS_DIR}
      RESULT_VARIABLE RESULT_CODE
    )

    IF(NOT RESULT_CODE EQUAL 0)
      MESSAGE( WARNING "Extracting the pre-compiled collage failed" )
    ENDIF()

  ENDIF()
ENDIF()

FILE(GLOB CollageEXE "${DOWNLOADED_APPS_DIR}/*")

IF(APPLE)
  INSTALL(FILES ${CollageEXE}
    DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
ELSE()
  INSTALL(FILES ${CollageEXE}
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
ENDIF()
