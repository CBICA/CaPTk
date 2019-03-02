# This CMakeLists file will add the source files in this directory for compilation

# Add all the files in here (or the full path of the directory) into cache so that they can be added as include files for applications
ADD_SUBDIRECTORY(src/IO)

PROJECT( CaPTk_CBICATK )

FILE( GLOB sources "${PROJECT_SOURCE_DIR}/src/*.*" )
SET( CBICA-TK_SOURCES ${sources} CACHE STRING "CBICA Toolkit sources" )
FILE( GLOB_RECURSE YAMLCPP_Headers "${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/src/*.h")
FILE( GLOB_RECURSE YAMLCPP_Sources "${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/src/*.cpp")
SET( YAMLCPP_Headers ${YAMLCPP_Headers} CACHE STRING "YAML-CPP headers" FORCE )
SET( YAMLCPP_Sources ${YAMLCPP_Sources} CACHE STRING "YAML-CPP sources" FORCE )
SOURCE_GROUP( YAMLCPP_Headers FILES ${YAMLCPP_Headers} )
SOURCE_GROUP( YAMLCPP_Sources FILES ${YAMLCPP_Sources} )

INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/src/
  ${PROJECT_SOURCE_DIR}/../thirdparty/
  ${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/include/ 
)
  
SET( CACHED_INCLUDE_DIRS 
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/src/
  ${PROJECT_SOURCE_DIR}/../thirdparty/
  ${PROJECT_SOURCE_DIR}/../thirdparty/yaml-cpp/include/ 
  CACHE STRING "All include directories"
)
  

SET( LIBNAME_CBICATK "${PROJECT_NAME}" CACHE STRING "CaPTk's common toolkit for use throughout the project" )
MARK_AS_ADVANCED( ${LIBNAME_CBICATK} )

IF(NOT (TARGET ${LIBNAME_CBICATK}))
  
  ADD_LIBRARY(
    ${LIBNAME_CBICATK}
    ${sources}
    ${YAMLCPP_Headers}
    ${YAMLCPP_Sources}
  )
  
  TARGET_LINK_LIBRARIES(
    ${LIBNAME_CBICATK}
    ${DEPENDENT_LIBS}
    CaPTk_IO
  )
  
ENDIF()


MARK_AS_ADVANCED(CBICA-TK_SOURCES)