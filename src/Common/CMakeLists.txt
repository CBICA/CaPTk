PROJECT(CaPTk_Common)

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY 
	${CaPTk_BINARY_DIR})
  
SET( ${PROJECT_NAME}_HDRS
  inc/CaPTkClassifierUtils.h
  inc/CaPTkDefines.h
  inc/CaPTkEnums.h
  inc/CaPTkGUIUtils.h
  inc/CaPTkUtils.h
)

SET( ${PROJECT_NAME}_SRCS
  src/CaPTkGUIUtils.cpp
)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/inc
                     ${CaPTk_SOURCE_DIR}/src/cbica_toolkit
                     ${CaPTk_SOURCE_DIR}/src/applications/common_includes
                     )
                     
ADD_LIBRARY(
  ${PROJECT_NAME}
  ${${PROJECT_NAME}_SRCS}
  ${${PROJECT_NAME}_HDRS}
)

