CMAKE_MINIMUM_REQUIRED(VERSION 3.7.2)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

SET( PROJECT_NAME Utilities )

CAPTK_ADD_PROJECT( ${PROJECT_NAME} 1.6.2 )

ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/ cbica_toolkit ) # takes care of find_itk

CAPTK_INITIAL_SETUP()

## dcmqi stuff starts

FIND_PACKAGE( DCMTK REQUIRED )

SET( CACHED_INCLUDE_DIRS
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/dcmqi/include/
  ${PROJECT_SOURCE_DIR}/dcmqi/include/dcmqi
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/json
)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/dcmqi/include/dcmqi)

set(HDRS
  ${INCLUDE_DIR}/preproc.h
  ${INCLUDE_DIR}/QIICRConstants.h
  ${INCLUDE_DIR}/QIICRUIDs.h
  ${INCLUDE_DIR}/ConverterBase.h
  ${INCLUDE_DIR}/Exceptions.h
  ${INCLUDE_DIR}/framesorter.h
  ${INCLUDE_DIR}/ImageSEGConverter.h
  ${INCLUDE_DIR}/ParaMapConverter
  ${INCLUDE_DIR}/Helper.h
  ${INCLUDE_DIR}/JSONMetaInformationHandlerBase.h
  ${INCLUDE_DIR}/JSONParametricMapMetaInformationHandler.h
  ${INCLUDE_DIR}/JSONSegmentationMetaInformationHandler.h
  ${INCLUDE_DIR}/SegmentAttributes.h
  ${INCLUDE_DIR}/TID1500Reader.h
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/json/json.h
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/json/json-forwards.h
  )

set(SRCS
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/ConverterBase.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/ImageSEGConverter.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/ParaMapConverter.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/Helper.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/JSONMetaInformationHandlerBase.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/JSONParametricMapMetaInformationHandler.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/JSONSegmentationMetaInformationHandler.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/SegmentAttributes.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/libsrc/TID1500Reader.cpp
  ${PROJECT_SOURCE_DIR}/dcmqi/jsoncpp/jsoncpp.cpp
  )


set(lib_name dcmqi)

add_library(${lib_name} STATIC
  ${HDRS}
  ${SRCS}
  ${ADDITIONAL_SRCS}
  )
  
target_link_libraries(${lib_name} PUBLIC
  ${DCMTK_LIBRARIES}
  ${ITK_LIBRARIES}
)
## dcmqi stuff ends

INCLUDE_DIRECTORIES( ${CACHED_INCLUDE_DIRS} )

CAPTK_ADD_EXECUTABLE( ${PROJECT_NAME} "${APPLICATION_SOURCES};${APPLICATION_DEPENDS}" ${LIBNAME_CBICATK} )

ADD_DEPENDENCIES( ${PROJECT_NAME} ${LIBNAME_CBICATK} ${lib_name} ) 