PROJECT(CaPTk_Common)

#SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY 
#	${CaPTk_BINARY_DIR})

SET( ${PROJECT_NAME}_HDRS
  ${PROJECT_SOURCE_DIR}/inc/CaPTkClassifierUtils.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkDefines.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkEnums.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkGUIUtils.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkUtils.h
  # ${PROJECT_SOURCE_DIR}/inc/ApplicationDownloadManager.h
  ${PROJECT_SOURCE_DIR}/inc/ApplicationPreferences.h
  # ${PROJECT_SOURCE_DIR}/inc/fAppDownloadDialog.h
  # ${PROJECT_SOURCE_DIR}/inc/ui_fAppDownloadDialog.h
  # ${PROJECT_SOURCE_DIR}/inc/QZipReader.h
  # ${PROJECT_SOURCE_DIR}/inc/QZipWriter.h
  # ${PROJECT_SOURCE_DIR}/inc/ThreadedExtraction.h
)

SET( ${PROJECT_NAME}_SRCS
  ${PROJECT_SOURCE_DIR}/src/CaPTkGUIUtils.cpp
  # ${PROJECT_SOURCE_DIR}/src/ApplicationDownloadManager.cpp
  ${PROJECT_SOURCE_DIR}/src/ApplicationPreferences.cpp
  # ${PROJECT_SOURCE_DIR}/src/fAppDownloadDialog.cpp
  # ${PROJECT_SOURCE_DIR}/src/QZip.cpp
  # ${PROJECT_SOURCE_DIR}/src/ThreadedExtraction.cpp
)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/inc
  ${CaPTk_SOURCE_DIR}/src/cbica_toolkit/src
  ${CaPTk_SOURCE_DIR}/src/applications/common_includes
)

ADD_LIBRARY(
  ${PROJECT_NAME}
  ${${PROJECT_NAME}_SRCS}
  ${${PROJECT_NAME}_HDRS}
)

TARGET_LINK_LIBRARIES(
  ${PROJECT_NAME}
  ${LIBNAME_CBICATK}
  ${DEPENDENT_LIBS}
)
 
SET( CACHED_INCLUDE_DIRS 
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/inc
  CACHE STRING "All include directories"
)
  