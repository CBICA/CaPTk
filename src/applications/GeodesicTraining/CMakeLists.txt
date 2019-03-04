CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

SET( PROJECT_NAME GeodesicTraining )

CAPTK_ADD_PROJECT( ${PROJECT_NAME} 1.0 )

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

# For C++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/ 
  ${PROJECT_SOURCE_DIR}/../../thirdparty/
  ${PROJECT_SOURCE_DIR}/../../thirdparty/yaml-cpp/include/ 
)

# --------- Library dependencies ---------

# For ITK
FIND_PACKAGE( ITK REQUIRED )
SET( ITK_NO_IO_FACTORY_REGISTER_MANAGER "OFF" )
INCLUDE( ${ITK_USE_FILE} )

# For OpenCV
FIND_PACKAGE( OpenCV REQUIRED )
#INCLUDE_DIRECTORIES( ${OpenCV_INCLUDE_DIRS} )

# For OpenMP
FIND_PACKAGE(OpenMP REQUIRED)
SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

FIND_PACKAGE( Qt5 COMPONENTS Core )

# --------- Source code dependencies ---------

FILE( GLOB_RECURSE CBICA_TOOLKIT_HEADERS
  "${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/*.h"
  "${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/IO/inc/*.h"
)
FILE( GLOB_RECURSE CBICA_TOOLKIT_SOURCES
  "${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/*.cpp"
  "${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/IO/src/*.cpp"
)
FILE( GLOB_RECURSE YAML_HEADERS
  "${PROJECT_SOURCE_DIR}/../../thirdparty/yaml-cpp/include/yaml-cpp/*.*"
)
FILE( GLOB_RECURSE YAML_SOURCES
  "${PROJECT_SOURCE_DIR}/../../thirdparty/yaml-cpp/src/*.*"
)


# --------- Setting the project files ---------

FILE( GLOB_RECURSE GEODESIC_TRAINING_HEADERS
  ${PROJECT_SOURCE_DIR}/include/*.h
)

FILE( GLOB_RECURSE GEODESIC_TRAINING_SOURCES
  ${PROJECT_SOURCE_DIR}/include/*.cpp
)

SET( GEODESIC_TRAINING_INCLUDE_DIRS 
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/../../thirdparty/thirdparty/yaml-cpp/include/yaml-cpp
  ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src
  ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/IO/inc
)


# --------- For the library ---------

ADD_LIBRARY( GeodesicTrainingCaPTkLib
  GeodesicTrainingCaPTkApp.cpp
  ${GEODESIC_TRAINING_SOURCES}
  ${CBICA_TOOLKIT_SOURCES}
  ${YAML_SOURCES}
)

TARGET_LINK_LIBRARIES( GeodesicTrainingCaPTkLib
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
)

qt5_use_modules( GeodesicTrainingCaPTkLib Core )

TARGET_INCLUDE_DIRECTORIES( GeodesicTrainingCaPTkLib PUBLIC
  ${GEODESIC_TRAINING_INCLUDE_DIRS}
)


# --------- For the executable ---------

SET( GEODESIC_TRAINING_EXE_FILES
  GeodesicTraining.cxx
  "${GEODESIC_TRAINING_HEADERS}"
  "${GEODESIC_TRAINING_SOURCES}"
  "${YAML_HEADERS}"
  "${YAML_SOURCES}"
  "${CBICA_TOOLKIT_HEADERS}"
  "${CBICA_TOOLKIT_SOURCES}"
)

SET( GEODESIC_TRAINING_EXE_LIBS
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
)

# CAPTK_ADD_EXECUTABLE( GeodesicTraining
#   ${GEODESIC_TRAINING_EXE_FILES}
#   ${GEODESIC_TRAINING_EXE_LIBS}
# )
