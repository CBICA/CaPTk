CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

SET( PROJECT_NAME GeodesicTrainingCaPTk )

# This cmake file is here because there were problems
# adding the GeodesicTraining cmake file directly

# Set project name 
PROJECT( ${PROJECT_NAME} )
SET( PROJECT_VERSION "0.12")
ADD_DEFINITIONS(-DPROJECT_VERSION="${PROJECT_VERSION}" )

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
 
# For ITK
FIND_PACKAGE( ITK REQUIRED )
SET(ITK_NO_IO_FACTORY_REGISTER_MANAGER "OFF")
INCLUDE( ${ITK_USE_FILE} )

# For OpenCV
FIND_PACKAGE( OpenCV REQUIRED )
#INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

# For Qt
find_package(Qt5 COMPONENTS Core Widgets REQUIRED)

# For OpenMP
FIND_PACKAGE(OpenMP REQUIRED)
SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

# For C++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

# Setting the files

FILE( GLOB_RECURSE AGD_HEADERS 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/AdaptiveGeodesicDistance/include/AdaptiveGeodesicDistance/*.h" 
)
FILE( GLOB_RECURSE AGD_SOURCES 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/AdaptiveGeodesicDistance/src/*.cpp" 
)
FILE( GLOB_RECURSE RF_HEADERS 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/RandomForestSuite/include/RandomForestSuite/*.h" 
)
FILE( GLOB_RECURSE RF_SOURCES 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/RandomForestSuite/src/*.cpp" 
)
FILE( GLOB_RECURSE SVM_HEADERS 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/include/SvmSuite/*.h" 
)
FILE( GLOB_RECURSE SVM_SOURCES 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/src/*.cpp" 
)
FILE( GLOB_RECURSE SVM_CONV_SOURCES 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/src/convertions/*.cpp" 
)
FILE( GLOB_RECURSE GTS_HEADERS 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/include/GeodesicTrainingSegmentation/*.h" 
)
FILE( GLOB_RECURSE GTS_SOURCES 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/src/*.cpp" 
)
FILE( GLOB GT_SOURCES 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/src/*.*"
)
#FILE( GLOB_RECURSE CBICA_TOOLKIT_HEADERS
#  "${PROJECT_SOURCE_DIR}/GeodesicTraining/src/cbica_toolkit/*.h"
#)
#FILE( GLOB_RECURSE CBICA_TOOLKIT_SOURCES
#  "${PROJECT_SOURCE_DIR}/GeodesicTraining/src/cbica_toolkit/*.cpp"
#)
FILE( GLOB_RECURSE CBICA_TOOLKIT_HEADERS
  "${CMAKE_SOURCE_DIR}/src/cbica_toolkit/src/*.h"
)
FILE( GLOB_RECURSE CBICA_TOOLKIT_SOURCES
  "${CMAKE_SOURCE_DIR}/src/cbica_toolkit/src/*cpp"
)
FILE( GLOB_RECURSE YAML_HEADERS
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/thirdparty/yaml-cpp/include/yaml-cpp/*.*"
)
FILE( GLOB_RECURSE YAML_SOURCES
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/thirdparty/yaml-cpp/src/*.*"
)

SET( AGD_HEADERS_DIR 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/AdaptiveGeodesicDistance/include/AdaptiveGeodesicDistance" 
)
SET( RF_HEADERS_DIR 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/RandomForestSuite/include/RandomForestSuite" 
)
SET( SVM_HEADERS_DIR 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/include/SvmSuite" 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/src/convertions"
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/SvmSuite/thirdparty/yaml-cpp/include"
)
SET( GTS_HEADERS_DIR 
  "${PROJECT_SOURCE_DIR}/GeodesicTraining/GeodesicTrainingSegmentation/include/GeodesicTrainingSegmentation" 
)
SET( CBICA_TOOLKIT_HEADERS_DIR
  "${CMAKE_SOURCE_DIR}/src/cbica_toolkit/src"
)

SET( GEODESIC_TRAINING_SOURCES
  ${AGD_SOURCES}
  ${RF_SOURCES}
  ${SVM_SOURCES}
  ${SVM_CONV_SOURCES}
  ${YAML_HEADERS}
  ${YAML_SOURCES}
  #${CBICA_TOOLKIT_SOURCES}
  ${GTS_SOURCES}
)

# These are needed to include the headers in CaPTk
SET( GEODESIC_TRAINING_HEADERS
  ${AGD_HEADERS}
  ${RF_HEADERS}
  ${SVM_HEADERS}
  ${GTS_HEADERS}
)

SET( GEODESIC_TRAINING_INCLUDE_DIRS 
  ${AGD_HEADERS_DIR}
  ${RF_HEADERS_DIR}
  ${SVM_HEADERS_DIR}
  ${GTS_HEADERS_DIR}
  ${PROJECT_SOURCE_DIR}
  ${CBICA_TOOLKIT_HEADERS_DIR}
  CACHE INTERNAL ""
)

# Creating the library
ADD_LIBRARY(${PROJECT_NAME} STATIC 
  GeodesicTrainingCaPTkApp.cpp
  ${GEODESIC_TRAINING_SOURCES}
)

# Linking the dependencies of the library
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
)

qt5_use_modules(${PROJECT_NAME} Core Widgets)

# Include folders for the library
TARGET_INCLUDE_DIRECTORIES( ${PROJECT_NAME} PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${GEODESIC_TRAINING_INCLUDE_DIRS}
)

SET( GEODESIC_TRAINING_EXE_FILES
  ${GT_SOURCES} 
  ${GEODESIC_TRAINING_HEADERS}
  ${GEODESIC_TRAINING_SOURCES}
  #${CBICA_TOOLKIT_HEADERS}
  #${CBICA_TOOLKIT_SOURCES}
  ${YAML_HEADERS}
  ${YAML_SOURCES}
)

SET( GEODESIC_TRAINING_EXE_LIBS
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
  ${LIBNAME_CBICATK}
)

# Add sources to executable
CAPTK_ADD_EXECUTABLE( GeodesicTraining
  "${GEODESIC_TRAINING_EXE_FILES}"
  "${GEODESIC_TRAINING_EXE_LIBS}"
)

# Include folders for the executable
TARGET_INCLUDE_DIRECTORIES( GeodesicTraining PRIVATE 
  ${GEODESIC_TRAINING_INCLUDE_DIRS}
)