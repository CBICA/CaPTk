# This CMakeLists file will add the source files in this directory for compilation

# Add all the files in here (or the full path of the directory) into cache so that they can be added as include files for applications

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

CAPTK_INITIAL_SETUP()

PROJECT( CaPTk_CBICATK )

FILE( GLOB_RECURSE sources_h "${PROJECT_SOURCE_DIR}/src/*.h" )
FILE( GLOB_RECURSE sources_hxx "${PROJECT_SOURCE_DIR}/src/*.hxx" )
FILE( GLOB_RECURSE sources_cpp "${PROJECT_SOURCE_DIR}/src/*.cpp" )

SET( CBICA-TK_SOURCES
  ${sources_h}
  ${sources_hxx}
  ${sources_cpp}
)

SET( CBICA-TK_SOURCES ${CBICA-TK_SOURCES} CACHE STRING "CBICA Toolkit sources" FORCE )
FILE( GLOB_RECURSE YAMLCPP_Headers "${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/includes/*.h")
FILE( GLOB_RECURSE YAMLCPP_Headers_MORE "${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/src/*.h")
FILE( GLOB_RECURSE YAMLCPP_Sources "${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/src/*.cpp")
SET( YAMLCPP_Headers ${YAMLCPP_Headers} ${YAMLCPP_Headers_MORE} CACHE STRING "YAML-CPP headers" FORCE )
SET( YAMLCPP_Sources ${YAMLCPP_Sources} CACHE STRING "YAML-CPP sources" FORCE )
SOURCE_GROUP( YAMLCPP_Headers FILES ${YAMLCPP_Headers} )
SOURCE_GROUP( YAMLCPP_Sources FILES ${YAMLCPP_Sources} )

#Find libraries
FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

SET( CACHED_INCLUDE_DIRS 
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/src/
  ${PROJECT_SOURCE_DIR}/src/IO/inc/
  ${PROJECT_SOURCE_DIR}/../thirdparty/
  ${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/include/ 
  CACHE STRING "All include directories" FORCE
)

INCLUDE_DIRECTORIES( 
  ${CACHED_INCLUDE_DIRS}
)
  
SET( LIBNAME_CBICATK "${PROJECT_NAME}" CACHE STRING "CaPTk's common toolkit for use throughout the project" )
MARK_AS_ADVANCED( ${LIBNAME_CBICATK} )

IF(NOT (TARGET ${LIBNAME_CBICATK}))
  
  ADD_LIBRARY(
    ${LIBNAME_CBICATK}
    ${CBICA-TK_SOURCES}
    ${YAMLCPP_Headers}
    ${YAMLCPP_Sources}
  )
  
  TARGET_LINK_LIBRARIES(
    ${LIBNAME_CBICATK}
    ${ITK_LIBRARIES}
    #CaPTk_IO
  )
  
ENDIF()

MARK_AS_ADVANCED(CBICA-TK_SOURCES)