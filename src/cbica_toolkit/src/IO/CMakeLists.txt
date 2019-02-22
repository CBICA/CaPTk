PROJECT(CaPTk_IO)

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY 
	${CaPTk_BINARY_DIR})
  
SET( ${PROJECT_NAME}_HDRS
  inc/DicomSeriesReader.h
  inc/DicomMetadataReader.h
)

SET( ${PROJECT_NAME}_SRCS
  src/DicomSeriesReader.cpp
  src/DicomMetadataReader.cpp
)
                     
INCLUDE_DIRECTORIES( 
${CMAKE_CURRENT_SOURCE_DIR}/inc
)

ADD_LIBRARY(
  ${PROJECT_NAME}
  ${${PROJECT_NAME}_SRCS}
  ${${PROJECT_NAME}_HDRS}
)

