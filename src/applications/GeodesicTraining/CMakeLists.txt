CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )


SET( PROJECT_NAME GeodesicTraining )

CAPTK_ADD_PROJECT( ${PROJECT_NAME} "1.0" )
ADD_DEFINITIONS( -DPROJECT_VERSION="1.0" ) # Needed

CAPTK_INITIAL_SETUP()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
  
INCLUDE_DIRECTORIES( 
  ${CACHED_INCLUDE_DIRS}
)

# --------- Library dependencies ---------

# For ITK
FIND_PACKAGE( ITK REQUIRED )
SET( ITK_NO_IO_FACTORY_REGISTER_MANAGER "OFF" )
INCLUDE( ${ITK_USE_FILE} )

# For OpenCV
FIND_PACKAGE( OpenCV REQUIRED )
#INCLUDE_DIRECTORIES( ${OpenCV_INCLUDE_DIRS} )

FIND_PACKAGE( Qt5 COMPONENTS Core )

# --------- Setting the project files ---------

FILE( GLOB_RECURSE GEODESIC_TRAINING_HEADERS
  ${PROJECT_SOURCE_DIR}/src/depends/*.h
)

FILE( GLOB_RECURSE GEODESIC_TRAINING_SOURCES
  ${PROJECT_SOURCE_DIR}/src/depends/*.cpp
)

ADD_LIBRARY( GeodesicTrainingComputeLib
  ${GEODESIC_TRAINING_HEADERS}
  ${GEODESIC_TRAINING_SOURCES}
)

TARGET_LINK_LIBRARIES( GeodesicTrainingComputeLib
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
  ${LIBNAME_CBICATK}
)

# --------- For the executable ---------

CAPTK_ADD_EXECUTABLE( ${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/GeodesicTraining.cxx "GeodesicTrainingComputeLib;${LIBNAME_CBICATK};${ITK_LIBRARIES}" )

ADD_DEPENDENCIES( ${PROJECT_NAME} ${LIBNAME_CBICATK} GeodesicTrainingComputeLib )

# --------- For the library ---------

ADD_LIBRARY( GeodesicTrainingCaPTkLib
  ${PROJECT_SOURCE_DIR}/src/GeodesicTrainingCaPTkApp.cpp
)

TARGET_LINK_LIBRARIES( GeodesicTrainingCaPTkLib
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
  ${LIBNAME_CBICATK}
  ${GeodesicTrainingComputeLib}
)

qt5_use_modules( GeodesicTrainingCaPTkLib Core )
