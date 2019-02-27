cmake_minimum_required(VERSION 3.0)
SET( PROJECT_NAME RandomForestSuite )
project (${PROJECT_NAME} VERSION 1.0 DESCRIPTION "Random Forest operations")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

# For OpenCV
FIND_PACKAGE( OpenCV REQUIRED )

# Setting paths for headers and sources
SET( Project_Include_Dir "${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}" )
SET( Project_Src_Dir "${PROJECT_SOURCE_DIR}/src" )

# Setting the sources
FILE(GLOB_RECURSE RF_SUITE_Sources "${Project_Src_Dir}/*.*")

# Including the paths for headers
include_directories(
    "${Project_Include_Dir}"
    "${OpenCV_INCLUDE_DIRS}"
)

# Creating the library
add_library(${PROJECT_NAME} STATIC 
    ${RF_SUITE_Sources}
)

# Linking the dependencies of the library
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
    ${OpenCV_LIBRARIES}
)
