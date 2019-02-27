CMAKE_MINIMUM_REQUIRED(VERSION 3.7.2)

SET( PROJECT_NAME FeatureExtraction )

# Set project name 
PROJECT( ${PROJECT_NAME} VERSION 1.6.1 LANGUAGES CXX )
#SET( PROJECT_VERSION "1.6.1")
ADD_DEFINITIONS(-DPROJECT_VERSION="${PROJECT_VERSION}" )

#MESSAGE( STATUS "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}" )
INCLUDE_DIRECTORIES( 
  ${PROJECT_SOURCE_DIR}/src/depends/ 
  ${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/ 
  ${PROJECT_SOURCE_DIR}/../../thirdparty/
  ${PROJECT_SOURCE_DIR}/../../thirdparty/yaml-cpp/include/ 
)

#Find libraries
FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

# For new feature, add the feature definition name here; it needs to correspond to the file in "${DATA_DIR}/features"
SET( FeatureDefinitions "" )
FILE( GLOB allFeatureDefinitionFiles "${PROJECT_SOURCE_DIR}/data/*.csv" )
FOREACH(featureDefinitionFile ${allFeatureDefinitionFiles})

  GET_FILENAME_COMPONENT( CURRENT_DEF_FILE ${featureDefinitionFile} NAME_WE )

  IF (NOT APPLE)
    INSTALL(
      FILES "${featureDefinitionFile}"
      DESTINATION data/features
    )
  ENDIF()

  SET( FeatureDefinitions "${FeatureDefinitions} ${CURRENT_DEF_FILE}" )
ENDFOREACH()

IF (APPLE)
  SET_SOURCE_FILES_PROPERTIES(${allFeatureDefinitionFiles} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/features)

  SET( BATCHMODE "${PROJECT_SOURCE_DIR}/data/features/batchMode/batch_featureExtraction.csv" )
  SET_SOURCE_FILES_PROPERTIES(${BATCHMODE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/share/featureExtractionBatch)

  SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${allFeatureDefinitionFiles} ${BATCHMODE})

ELSE()
  INSTALL( FILES "${PROJECT_SOURCE_DIR}/data/batchMode/batch_featureExtraction.csv" DESTINATION share/featureExtractionBatch )
ENDIF()

SET( FeatureDefinitions ${FeatureDefinitions} CACHE STRING "Available Feature Definitions" FORCE)
ADD_DEFINITIONS(-DFeatureDefinitions="${FeatureDefinitions}")

if(APPLE)
#  set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/6.0.1/bin/clang")
#  set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/6.0.1/bin/clang++")
#  set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/6.0.1/lib")
#  set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/6.0.1/include")
#
#  OPTION (USE_OpenMP "Use OpenMP to enable <omp.h>" ON)
#
#  # Find OpenMP
#  if(USE_OpenMP)
#      if(CMAKE_C_COMPILER_ID MATCHES "Clang")
#          set(OpenMP_C "${CMAKE_C_COMPILER}")
#          set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
#          set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
#          set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
#          set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
#          set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
#      endif()
#      if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#          set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
#          set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
#          set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
#          set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
#          set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
#          set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
#      endif()
#  endif()
#
#  if(USE_OpenMP)
#      find_package(OpenMP REQUIRED)
#  endif(USE_OpenMP)
#
#  if (OPENMP_FOUND)
#    include_directories("${OPENMP_INCLUDES}")
#    link_directories("${OPENMP_LIBRARIES}")
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#  endif(OPENMP_FOUND)
else()
  FIND_PACKAGE(OpenMP REQUIRED)
  SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 
SET_PROPERTY( GLOBAL PROPERTY USE_FOLDERS ON )

FILE( GLOB sources "${PROJECT_SOURCE_DIR}/src/depends/*.*" )
FILE( GLOB cbicatk "${PROJECT_SOURCE_DIR}/../../cbica_toolkit/src/*.*" )
FILE( GLOB_RECURSE YAMLCPP_Headers "${PROJECT_SOURCE_DIR}/../../thirdparty/yaml-cpp/src/*.h")
FILE( GLOB_RECURSE YAMLCPP_Sources "${PROJECT_SOURCE_DIR}/../../thirdparty/yaml-cpp/src/*.cpp")
SET( YAMLCPP_Headers ${YAMLCPP_Headers} CACHE STRING "YAML-CPP headers" FORCE )
SET( YAMLCPP_Sources ${YAMLCPP_Sources} CACHE STRING "YAML-CPP sources" FORCE )
SOURCE_GROUP( YAMLCPP_Headers FILES ${YAMLCPP_Headers} )
SOURCE_GROUP( YAMLCPP_Sources FILES ${YAMLCPP_Sources} )

SET( LIBNAME_CBICATK "CaPTk_CBICATK" )
IF(NOT (TARGET ${LIBNAME_CBICATK}))
  
  ADD_LIBRARY(
    ${LIBNAME_CBICATK}
    ${cbicatk}
    ${YAMLCPP_Sources}
  )
  
  TARGET_LINK_LIBRARIES(
    ${LIBNAME_CBICATK}
    ${ITK_LIBRARIES}
    ${OpenCV_LIBRARIES}
  )
  
ENDIF()

CAPTK_ADD_EXECUTABLE(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/FeatureExtraction.cxx "${LIBNAME_CBICATK}")

ADD_DEPENDENCIES( ${PROJECT_NAME} ${LIBNAME_CBICATK} )

# Link the libraries to be used
TARGET_LINK_LIBRARIES(
  ${PROJECT_NAME}
  ${LIBNAME_CBICATK}
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
)

OPTION( BUILD_REFERENCES "Build reference executables from CBICA pipeline" OFF )

IF( BUILD_REFERENCES )

  SET_PROPERTY( GLOBAL PROPERTY USE_FOLDERS ON )
  SET( StandAloneCLIAppsFolder "References")
  # find stand-alone CMD applications
  FILE( GLOB STANDALONE_APPS
    "${PROJECT_SOURCE_DIR}/src/reference/*.cxx"
  )

  INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/src/reference/ )
    
  FOREACH(application ${STANDALONE_APPS})

    GET_FILENAME_COMPONENT( STANDALONE_APP_NAME ${application} NAME_WE )

    ADD_EXECUTABLE( 
      ${STANDALONE_APP_NAME}
      ${PROJECT_SOURCE_DIR}/src/reference/${STANDALONE_APP_NAME}.cxx
    )
    
    TARGET_LINK_LIBRARIES( 
      ${STANDALONE_APP_NAME}
      ${ITK_LIBRARIES}
      ${OpenCV_LIBRARIES}
    )
    
    SET_TARGET_PROPERTIES( ${STANDALONE_APP_NAME} PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )
    
  ENDFOREACH()
  
ENDIF()

#INSTALL( FILES ${PROJECT_SOURCE_DIR}/data/1_params_default.csv DESTINATION data/features )
#INSTALL( FILES ${PROJECT_SOURCE_DIR}/data/2_params_default_lattice.csv DESTINATION data/features )
#INSTALL( FILES ${PROJECT_SOURCE_DIR}/data/batchMode/batch_featureExtraction.csv DESTINATION data/features/batchMode )

# For new feature, add the feature definition name here; it needs to correspond to the file in "${DATA_DIR}/features"
SET( FeatureDefinitions "" )
FILE( GLOB allFeatureDefinitionFiles "${PROJECT_SOURCE_DIR}/data/*.csv" )
FOREACH(featureDefinitionFile ${allFeatureDefinitionFiles})

  GET_FILENAME_COMPONENT( CURRENT_DEF_FILE ${featureDefinitionFile} NAME_WE )

  IF( APPLE )
    SET_SOURCE_FILES_PROPERTIES(${featureDefinitionFile} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/features)
  ELSE()
    INSTALL(
      FILES "${featureDefinitionFile}"
      DESTINATION data/features
    )
  ENDIF()

  SET( FeatureDefinitions "${FeatureDefinitions} ${CURRENT_DEF_FILE}" )
ENDFOREACH()

IF (APPLE)
  SET_SOURCE_FILES_PROPERTIES(${allFeatureDefinitionFiles} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/features)

  SET( BATCHMODE "${PROJECT_SOURCE_DIR}/data/batchMode/batch_featureExtraction.csv" )
  SET_SOURCE_FILES_PROPERTIES(${BATCHMODE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/share/featureExtractionBatch)

  SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${allFeatureDefinitionFiles} ${BATCHMODE})

ELSE()
  INSTALL( FILES "${PROJECT_SOURCE_DIR}/data/batchMode/batch_featureExtraction.csv" DESTINATION share/featureExtractionBatch )
ENDIF()

SET( FeatureDefinitions ${FeatureDefinitions} CACHE STRING "Available Feature Definitions" FORCE)
ADD_DEFINITIONS(-DFeatureDefinitions="${FeatureDefinitions}")
