PROJECT(CaPTk_Common)

#SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY 
#	${CaPTk_BINARY_DIR})

SET(CMAKE_AUTOMOC ON)

SET( ${PROJECT_NAME}_HDRS
  ${PROJECT_SOURCE_DIR}/inc/CaPTkClassifierUtils.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkDefines.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkEnums.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkGUIUtils.h
  ${PROJECT_SOURCE_DIR}/inc/CaPTkUtils.h
  ${PROJECT_SOURCE_DIR}/inc/DownloadManager.h
)

SET( ${PROJECT_NAME}_SRCS
  ${PROJECT_SOURCE_DIR}/src/CaPTkGUIUtils.cpp
  ${PROJECT_SOURCE_DIR}/src/DownloadManager.cpp
)

#Message(STATUS "captk src dir:" ${CaPTk_SOURCE_DIR})

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/inc
  ${CaPTk_SOURCE_DIR}/src/cbica_toolkit/src
  ${PROJECT_SOURCE_DIR}/src/applications/common_includes
)

ADD_LIBRARY(
  ${PROJECT_NAME}
  ${${PROJECT_NAME}_SRCS}
  ${${PROJECT_NAME}_HDRS}
)

TARGET_LINK_LIBRARIES(
  ${PROJECT_NAME}
    Qt5::Core
    Qt5::Gui
    Qt5::Widgets
    Qt5::Concurrent
    Qt5::Network
 # ${LIBNAME_CBICATK}
 # ${DEPENDENT_LIBS}
)
 
SET( CACHED_INCLUDE_DIRS 
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/inc
  CACHE STRING "All include directories"
)
  
