cmake_minimum_required(VERSION 3.0)
SET( PROJECT_NAME AdaptiveGeodesicDistance )
project (${PROJECT_NAME} VERSION 1.0 DESCRIPTION "Adaptive Geodesic Distance Algorithm")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

# Setting paths for headers and sources
SET( Project_Include_Dir "${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}" )
SET( Project_Src_Dir "${PROJECT_SOURCE_DIR}/src" )

# Setting the sources
FILE(GLOB_RECURSE AGD_Sources "${Project_Src_Dir}/*.*")

# Including the path for headers
include_directories("${Project_Include_Dir}")

# For ITK
FIND_PACKAGE( ITK REQUIRED )
SET(ITK_NO_IO_FACTORY_REGISTER_MANAGER "OFF")
INCLUDE( ${ITK_USE_FILE} )

# Creating the library
add_library(${PROJECT_NAME} STATIC 
	${AGD_Sources}
)

# Linking the dependencies of the library
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
    ${ITK_LIBRARIES}
)
