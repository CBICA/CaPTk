# This CMakeLists file will add the source files in this directory for compilation

# Add all the files in here (or the full path of the directory) into cache so that they can be added as include files for applications

SET( CBICA-TK_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaUtilities.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaUtilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaCmdParser.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaCmdParser.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaLogging.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaLogging.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaITKImageInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaITKImageInfo.cpp
  
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaITKUtilities.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaITKSafeImageIO.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaDTIProcessingManager.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cbicaDTIProcessingManager.cpp
  CACHE STRING "CBICA Toolkit sources"
)

MARK_AS_ADVANCED(CBICA-TK_SOURCES)