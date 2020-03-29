CMAKE_MINIMUM_REQUIRED(VERSION 3.7.2)

# set project name. The exe name is taken from this
SET( PROJECT_NAME CaPTk )
PROJECT( CaPTk )
SET( ${PROJECT_NAME} )
SET( ${PROJECT_NAME}_Variant "Full" ) # the particular variant of CaPTk (Full/Neuro/Torso/Viewer)

# WARNING: captk-pkg (the linux packager) directly tries to pattern match these lines (excluding the actual number/patch text)
# to read the version number.
# So small changes here to formatting or naming convention *may potentially break the linux packager*!
# STRING(TIMESTAMP TODAY "%Y%m%d")

SET( PROJECT_VERSION_MAJOR 1 )
SET( PROJECT_VERSION_MINOR 7 )
SET( PROJECT_VERSION_PATCH 7.nonRelease )
SET( PROJECT_VERSION_TWEAK )

# check for the string "nonRelease" in the PROJECT_VERSION_PATCH variable
# if it is found, add timestamp and git hash of commit to the version
STRING(FIND ${PROJECT_VERSION_PATCH} "nonRelease" NONRELEASE_FOUND)
IF( ${NONRELEASE_FOUND} GREATER 0 )
  # get the timestamp
  STRING(TIMESTAMP TODAY "%Y%m%d")

  # get git hash
  IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    FIND_PACKAGE(Git)
    IF(GIT_FOUND)
      EXECUTE_PROCESS(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE "GIT_HASH"
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      #MESSAGE( STATUS "Git version: ${GIT_HASH}" )
    ELSE(GIT_FOUND)
      SET(GIT_HASH 0)
    ENDIF(GIT_FOUND)
  ENDIF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  
  SET( PROJECT_VERSION_PATCH "${PROJECT_VERSION_PATCH}.${TODAY}.${GIT_HASH}" )
  #MESSAGE( STATUS "PROJECT_VERSION_PATCH: ${PROJECT_VERSION_PATCH}")

ENDIF()

# Configure CCache if available
FIND_PROGRAM(CCACHE_FOUND ccache)
IF(CCACHE_FOUND)
  SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
ENDIF(CCACHE_FOUND)

FIND_PROGRAM(CCacheExectuable ccache)
IF( CCacheExectuable )
  # only used with >=cmake-3.4
  SET( CMAKE_C_COMPILER_LAUNCHER   "${CCacheExectuable}" )
  SET( CMAKE_CXX_COMPILER_LAUNCHER "${CCacheExectuable}" )
  #MESSAGE( STATUS "Speeding up build with ccache" )
ENDIF()

IF(WIN32)
  # enable parallel builds for windows
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
ENDIF()

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules/)
INCLUDE( CaPTk_macros )

#SET (SVN_REVISION "${SVN_REVISION}")

OPTION( CAPTK_SUPERBUILD_FORCE "Force CaPTk Superbuild regardless of whether dependencies are found or not" OFF )

# claire option is not present for Windows, unfortunately
IF( UNIX )
  SET( BUILD_CLAIRE OFF CACHE BOOL "Do you want to build CLAIRE (this changes the Superbuild)?" )
ENDIF()

SET( QT_VERSION_ACTUAL "5.12.1" )
if(UNIX)
  if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    if(EXISTS "/etc/issue")
      file(READ "/etc/issue" LINUX_ISSUE)
      if(LINUX_ISSUE MATCHES "CentOS")
        SET( QT_VERSION_ACTUAL "5.11.1" )
      endif()
    endif()
  endif()
endif()

SET( QT_VERSION ${QT_VERSION_ACTUAL} CACHE STRING "Superbuild Qt version to build CaPTk against" FORCE )
SET( OPENCV_VERSION "3.4.7" CACHE STRING "Superbuild OpenCV version to build CaPTk against" FORCE )

SET(CAPTK_BUILD_DIR "${PROJECT_BINARY_DIR}" CACHE STRING "CaPTk build directory" FORCE )
ADD_DEFINITIONS(-DCAPTK_BUILD_DIR="${PROJECT_BINARY_DIR}")

CAPTK_INITIAL_SETUP()

IF( CAPTK_SUPERBUILD_FORCE )

  PROJECT(${PROJECT_NAME}_Depends)
  INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/Superbuild.cmake )
  return()

ELSE()
  FIND_PACKAGE(ITK)
  
  IF(NOT EXISTS ${ITK_USE_FILE})
    #MESSAGE( AUTHOR_WARNING "=== [DEBUG] PROJECT_BIN_DIR=${PROJECT_BINARY_DIR}" )
    IF( NOT EXISTS ${PROJECT_BINARY_DIR}/ITK-build/ITKConfig.cmake )

      OPTION( CAPTK_SUPERBUILD "Do a Superbuild of the CaPTk dependencies or Specify ITK install directory with the defined components" ON )
      IF( CAPTK_SUPERBUILD )
        PROJECT(${PROJECT_NAME}_Depends)
        INCLUDE( ${PROJECT_SOURCE_DIR}/cmake_modules/Superbuild.cmake )
        return()
      ENDIF()

    ENDIF()
  ENDIF()
  
  IF( EXISTS ${PROJECT_BINARY_DIR}/ITK-build/ITKConfig.cmake )
    SET( ITK_DIR "${PROJECT_BINARY_DIR}/ITK-build/" CACHE STRING "ITK_DIR Path for use in other builds" FORCE )
    SET( DCMTK_DIR "${PROJECT_BINARY_DIR}/DCMTK-build/" CACHE STRING "ITK_DIR Path for use in other builds" FORCE )

    IF (APPLE) # THIS QT FOLDER STRUCTURE WAS NOT REQUIRED FOR 5.11.2 
      SET( ENV{PATH} "${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/clang_64/lib/cmake/Qt5/:${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/clang_64/bin:$ENV{PATH}" CACHE PATH "" FORCE)
      SET( ENV{CMAKE_PREFIX_PATH} "${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/clang_64/lib/cmake/Qt5/:${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/clang_64/bin:${PROJECT_BINARY_DIR}/ITK-build/:${PROJECT_BINARY_DIR}/DCMTK-build/:${PROJECT_BINARY_DIR}/OpenCV-build/:${PROJECT_BINARY_DIR}/VTK-build:$ENV{CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE )
      SET( ENV{CMAKE_PROGRAM_PATH} "${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/clang_64/lib/cmake/Qt5/:${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/clang_64/bin:$ENV{CMAKE_PROGRAM_PATH}" CACHE PATH "" FORCE )
    ELSE()
      SET( ENV{CMAKE_PREFIX_PATH} "${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/lib/cmake/Qt5/:${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/bin:${PROJECT_BINARY_DIR}/ITK-build/:${PROJECT_BINARY_DIR}/DCMTK-build/:${PROJECT_BINARY_DIR}/OpenCV-build/:${PROJECT_BINARY_DIR}/VTK-build:$ENV{CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE )
      SET( ENV{LD_LIBRARY_PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/:$ENV{LD_LIBRARY_PATH}" CACHE PATH "" FORCE )
      SET( ENV{CMAKE_PROGRAM_PATH} "${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/lib/cmake/Qt5/:${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/bin:$ENV{CMAKE_PROGRAM_PATH}" CACHE PATH "" FORCE )
      SET( ENV{PATH} "${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/lib/cmake/Qt5/:${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/bin:$ENV{PATH}" CACHE PATH "" FORCE)
    ENDIF()

    #MESSAGE( AUTHOR_WARNING "=== [DEBUG] ENV{CMAKE_PREFIX_PATH}=$ENV{CMAKE_PREFIX_PATH}" )

  ENDIF()
ENDIF()

SET( PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}" CACHE STRING "The project version" FORCE )

MESSAGE( AUTHOR_WARNING "Building CaPTk version ${PROJECT_VERSION}")

IF (APPLE)
  SET (EXE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}")
ELSE()
  SET( EXE_NAME "${PROJECT_NAME}" )
ENDIF()

#IF( NOT( ${SVN_REVISION} STREQUAL "" ) )
#  SET( PROJECT_VERSION "${PROJECT_VERSION}_${SVN_REVISION}")
#ENDIF()
IF(${PROJECT_NAME}_Variant MATCHES "Full" )
  SET( PROJECT_NAME_EXTENDED "Cancer Imaging Phenomics Toolkit (${PROJECT_NAME})")
ELSE()
  SET( PROJECT_NAME_EXTENDED "Cancer Imaging Phenomics Toolkit (${PROJECT_NAME}_${${PROJECT_NAME}_Variant})" )
ENDIF()
SET( PROJECT_VENDOR "CBICA - UPenn" )
SET( LICENSE_FILE "${PROJECT_SOURCE_DIR}/LICENSE" )
SET( README_FILE  "${PROJECT_SOURCE_DIR}/README.txt" )
SET( PROJECT_DOMAIN "https://www.med.upenn.edu/cbica/captk/" )
SET( PROJECT_CONTACT "software@cbica.upenn.edu" )
SET( COPYRIGHT "Copyright (c) 2019 CBICA. All rights reserved.")
SET( IDENTIFIER "")
ADD_DEFINITIONS(-DPROJECT_NAME="${PROJECT_NAME}" )
ADD_DEFINITIONS(-DPROJECT_VERSION="${PROJECT_VERSION}" )
ADD_DEFINITIONS(-DEXE_NAME="${EXE_NAME}")
ADD_DEFINITIONS(-DCaPTk_Variant="${CaPTk_Variant}")
ADD_DEFINITIONS(-DPROJECT_NAME_EXTENDED="${PROJECT_NAME_EXTENDED}")
ADD_DEFINITIONS(-DCONTACT="${PROJECT_CONTACT}")
ADD_DEFINITIONS(-DPROJECT_CONTACT="${PROJECT_CONTACT}")
ADD_DEFINITIONS(-DPROJECT_DOMAIN="${PROJECT_DOMAIN}")
ADD_DEFINITIONS(-DPROJECT_LICENSE="${LICENSE_FILE}")
ADD_DEFINITIONS(-DLICENSE="${LICENSE_FILE}")
ADD_DEFINITIONS(-DPROJECT_README="${README_FILE}")
ADD_DEFINITIONS(-DREADME="${README_FILE}")
ADD_DEFINITIONS(-DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}/")
ADD_DEFINITIONS(-DPROJECT_VARIANT="${${PROJECT_NAME}_Variant}")

### Building the documentation
## The following 4 lines are for generation of public documentation
SET( BUILD_DOCUMENTATION ON CACHE BOOL "Do you want to build Documentation (needs Doxygen)?" )
SET( BUILD_DOCUMENTATION_ONLY OFF CACHE BOOL "This doesn't build anything else in the project" )
SET( CAPTK_UPDATE_DOCS OFF CACHE BOOL "Copy HTML files to $PROJECT_SOURCE_DIR/docs for git commit" )
MARK_AS_ADVANCED(BUILD_DOCUMENTATION)
MARK_AS_ADVANCED(BUILD_DOCUMENTATION_ONLY)
IF( BUILD_DOCUMENTATION_ONLY )
  SET( BUILD_DOCUMENTATION ON )
  SET( CAPTK_UPDATE_DOCS ON ) #in this case, the HTML files are updated all the time 
ENDIF()
IF( BUILD_DOCUMENTATION )
  ADD_SUBDIRECTORY( docs_sources )
ENDIF ()

IF( NOT BUILD_DOCUMENTATION_ONLY )
  # Provide build configuration option
  SET( CAPTK_BUILD_OPTION "PACKAGE_PROJECT" CACHE STRING "Toggle between packaging mode (optimized code) or developer mode (provides debugging info)" )

  SET_PROPERTY( CACHE CAPTK_BUILD_OPTION PROPERTY STRINGS PACKAGE_PROJECT DEVELOPER_MODE )
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

  IF( CAPTK_BUILD_OPTION STREQUAL PACKAGE_PROJECT)

    SET( CAPTK_PACKAGE_PROJECT ON )
    SET( CAPTK_UPDATE_DOCS ON CACHE BOOL "Copy HTML files to $PROJECT_SOURCE_DIR/docs for git commit" )
    SET( CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Default configuration types" FORCE )
    ADD_DEFINITIONS(-DCAPTK_PACKAGE_PROJECT)

    # add a few compiler options for GCC and MSVC
    IF( CMAKE_COMPILER_IS_GNUCXX )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -O3" )
    ENDIF()

    if(APPLE)

      SET (CAPTK_MACOSX_BUNDLE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}")
      SET (CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${PROJECT}" )
      
    endif()

    IF(MSVC)
      ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
    ENDIF(MSVC)

  ELSE()

    SET( CAPTK_PACKAGE_PROJECT OFF )
    SET( CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "Default configuration types" FORCE )

  ENDIF()

  OPTION( CAPTK_DISABLE_CWL "Disable CWL file generation - useful for debug builds or devops" ON )
  MARK_AS_ADVANCED(CAPTK_PACKAGE_PROJECT)
  MARK_AS_ADVANCED(CAPTK_BUILD_OPTION)
  MARK_AS_ADVANCED(CAPTK_DISABLE_CWL)

  # MESSAGE( AUTHOR_WARNING "=== [DEBUG] Building CaPTk version=${PROJECT_VERSION}" )
  
  SET( CMAKE_COLOR_MAKEFILE ON )
  # option to compile only the viewer
  SET( CAPTK_BUILD_CONSOLE_ONLY OFF CACHE BOOL "Do you want only the interactive capabilities of ${PROJECT_NAME}? Enabling this will remove all algorithms from compilation." )

  IF (CAPTK_BUILD_CONSOLE_ONLY)
    ADD_DEFINITIONS(-DCAPTK_BUILD_CONSOLE_ONLY)
  ENDIF()

  MARK_AS_ADVANCED(CAPTK_BUILD_CONSOLE_ONLY)

  SET( CAPTK_CLI_MODE OFF CACHE BOOL "Builds only CLI apps of CaPTk." )
  ADD_DEFINITIONS(-DCAPTK_CLI_MODE)

  # find dependencies START

  SET( THIRDPARTY_LIB_PATH ${PROJECT_SOURCE_DIR}/src/thirdparty )
  INCLUDE_DIRECTORIES( ${THIRDPARTY_LIB_PATH} )

  ### third party library processing
  # find and include all third party libraries
  SUBDIRLIST(SUBDIRECTORIES ${THIRDPARTY_LIB_PATH})
  
  FIND_PACKAGE( ITK REQUIRED )
  #FIND_PACKAGE( ITK REQUIRED )
  INCLUDE("${ITK_USE_FILE}")

  FIND_PACKAGE(OpenCV REQUIRED)
  INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

  FOREACH(thirdParty_include ${SUBDIRECTORIES})
    INCLUDE_DIRECTORIES( ${THIRDPARTY_LIB_PATH}/${thirdParty_include} )
    IF( EXISTS "${THIRDPARTY_LIB_PATH}/${thirdParty_include}/include" )
      INCLUDE_DIRECTORIES( ${THIRDPARTY_LIB_PATH}/${thirdParty_include}/include )
    ENDIF()
    #SUBDIRLIST(SUBDIRECTORIES_INT ${THIRDPARTY_LIB_PATH})
    #FOREACH(thirdParty_int_include ${THIRDPARTY_LIB_PATH}/${thirdParty_include})
    #  INCLUDE_DIRECTORIES( ${THIRDPARTY_LIB_PATH}/${thirdParty_include}/${thirdParty_int_include} )
    #ENDFOREACH()
  ENDFOREACH()

  # for YAML-CPP
  FILE(GLOB_RECURSE YAMLCPP_Headers "${PROJECT_SOURCE_DIR}/src/thirdparty/yaml-cpp/include/yaml-cpp/*.h")
  FILE(GLOB_RECURSE YAMLCPP_Sources "${PROJECT_SOURCE_DIR}/src/thirdparty/yaml-cpp/src/*.cpp")
  SET( YAMLCPP_Headers ${YAMLCPP_Headers} CACHE STRING "YAML-CPP headers" FORCE )
  SET( YAMLCPP_Sources ${YAMLCPP_Sources} CACHE STRING "YAML-CPP sources" FORCE )
  ### third party library processing

  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE("${VTK_USE_FILE}")
  IF(VTK_VERSION VERSION_LESS 6.0.0)
    SET(VTK_LIBRARIES
      vtkCommon
      vtkRendering
      vtkIO
      vtkFiltering
      vtkGraphics
      vtkWidgets
      vtkImaging
      vtkHybrid
      QVTK
    )
  ENDIF(VTK_VERSION VERSION_LESS 6.0.0)
  IF(VTK_VERSION VERSION_LESS 5.6.0)
    SET(VTK_LIBRARIES
      ${VTK_LIBRARIES}
      #vtkQVTK
    )
  ENDIF(VTK_VERSION VERSION_LESS 5.6.0)


  FIND_PACKAGE(Qt5 COMPONENTS Core Gui Svg Widgets WebView WebEngine WebEngineCore WebEngineWidgets Concurrent Test Sql REQUIRED)

  # For some reason IF(LINUX) doesn't work properly on linux
  # This is necessary for the whole project to compile on linux
  IF(WIN32)

  ELSEIF(APPLE)

  ELSE()
    FIND_PACKAGE(Qt5 COMPONENTS X11Extras REQUIRED)
  ENDIF()

  LINK_DIRECTORIES(${QT_LIBRARY_DIR})
  #INCLUDE(${QT_USE_FILE})
  SET(QT_USE_QTGUI TRUE)
  SET(QT_USE_QTMAIN TRUE)
  SET(QT_USE_QTGUI TRUE)
  ADD_DEFINITIONS(${QT_DEFINITIONS})
  ADD_DEFINITIONS(-DQT_PLUGIN)
  ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)
  LINK_DIRECTORIES(${QT_LIBRARY_DIR})

  # link dependent libraries
  SET( DEPENDENT_LIBS
    Qt5::Core
    Qt5::Gui
    Qt5::Svg
    Qt5::Widgets
    Qt5::WebView
    Qt5::WebEngine
    Qt5::WebEngineCore
    Qt5::WebEngineWidgets
    Qt5::Concurrent
    Qt5::Test
    Qt5::Sql
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
    ${OpenCV_LIBRARIES}
    ${OpenMP_LIBRARIES}
  )

  IF(UNIX AND NOT APPLE)
    SET( DEPENDENT_LIBS 
      ${DEPENDENT_LIBS}
      Qt5::X11Extras
    )
  ENDIF()
  
  SET( DEPENDENT_LIBS ${DEPENDENT_LIBS} CACHE STRING "All dependent libraries" )

  # find dependencies END

  SET( DATA_DIR "${PROJECT_SOURCE_DIR}/data" )

  # define build options for applications and add the respective definitions
  SET( BUILD_TRAINING ON CACHE BOOL "Do you want to build Training Module?" )
  SET( BUILD_EGFRvIII ON CACHE BOOL "Do you want to build EGFRvIII Estimation?" )
  SET( BUILD_EGFRvIIISVM ON CACHE BOOL "Do you want to build EGFRvIII SVM Estimation?" )
  SET( BUILD_RECURRENCE ON CACHE BOOL "Do you want to build Recurrence Estimation?" )
  SET( BUILD_PSEUDOPROGRESSION ON CACHE BOOL "Do you want to build Pseudoprogression Estimation?" )
  SET( BUILD_SURVIVAL ON CACHE BOOL "Do you want to build Survival Estimation?" )
  SET( BUILD_ISUBTYPE OFF CACHE BOOL "Do you want to build Imaging Subtype Estimation?" )
  SET( BUILD_MSUBTYPE ON CACHE BOOL "Do you want to build Molecular Subtype Estimation?" )
  SET( BUILD_ATLAS ON CACHE BOOL "Do you want to build Population Atlas?" )
  SET( BUILD_ALIGNMENT ON CACHE BOOL "Do you want to build Perfusion Alignment?" )
  SET( BUILD_PCA ON CACHE BOOL "Do you want to build Perfusion PCA?" )
  SET( BUILD_GEODESIC ON CACHE BOOL "Do you want to build Geodesic Segmentation?" )
  SET( BUILD_GEODESICTRAINING ON CACHE BOOL "Do you want to build Geodesic Training Segmentation?" )
  SET( BUILD_SBRT ON CACHE BOOL "Do you want to build SBRT Lung?" )
  SET( BUILD_WHITESTRIPE ON CACHE BOOL "Do you want to build WhiteStripe Normalization?" )
  SET( BUILD_ITKSNAP ON CACHE BOOL "Do you want to connect ITK-SNAP to CaPTk Console?" )
  SET( BUILD_LIBRASINGLE ON CACHE BOOL "Do you want to connect LIBRA to CaPTk Console (needs MATLAB runtime)?" )
  SET( BUILD_LIBRABATCH ON CACHE BOOL "Do you want to connect LIBRA to CaPTk Console (needs MATLAB runtime)?" )
  SET( BUILD_CONFETTI ON CACHE BOOL "Do you want to connect Confetti to CaPTk Console (needs Python 2.7 and other dependencies)?" )
  SET( BUILD_FETALBRAIN OFF CACHE BOOL "Do you want to build FetalBrain?" )
  SET( BUILD_DEEPMEDIC ON CACHE BOOL "Do you want to package DeepMedic?" )

  IF(BUILD_EGFRvIII)
    ADD_DEFINITIONS(-DBUILD_EGFRvIII)
  ENDIF()

  IF(BUILD_EGFRvIIISVM)
    ADD_DEFINITIONS(-DBUILD_EGFRvIIISVM)
  ENDIF()
  
  IF(BUILD_ALIGNMENT)
    ADD_DEFINITIONS(-DBUILD_ALIGNMENT)
  ENDIF()

  IF(BUILD_PCA)
    ADD_DEFINITIONS(-DBUILD_PCA)
  ENDIF()
  
  IF(BUILD_RECURRENCE)
    ADD_DEFINITIONS(-DBUILD_RECURRENCE)
  ENDIF()

  IF(BUILD_PSEUDOPROGRESSION)
    ADD_DEFINITIONS(-DBUILD_PSEUDOPROGRESSION)
  ENDIF()

  IF(BUILD_ISUBTYPE)
    ADD_DEFINITIONS(-DBUILD_ISUBTYPE)
  ENDIF()

  IF(BUILD_MSUBTYPE)
    ADD_DEFINITIONS(-DBUILD_MSUBTYPE)
  ENDIF()

  IF(BUILD_ATLAS)
    ADD_DEFINITIONS(-DBUILD_ATLAS)
  ENDIF()

  IF(BUILD_RECURRENCE)
    ADD_DEFINITIONS(-DBUILD_RECURRENCE)
  ENDIF()

  IF(BUILD_SURVIVAL)
    ADD_DEFINITIONS(-DBUILD_SURVIVAL)
  ENDIF()

  IF(BUILD_GEODESIC)
    ADD_DEFINITIONS(-DBUILD_GEODESIC)
  ENDIF()

  IF(BUILD_SBRT)
    ADD_DEFINITIONS(-DBUILD_SBRT)
  ENDIF()

  IF(BUILD_WHITESTRIPE)
    ADD_DEFINITIONS(-DBUILD_WHITESTRIPE)
  ENDIF()

  IF(BUILD_ITKSNAP)
    ADD_DEFINITIONS(-DBUILD_ITKSNAP)
  ENDIF()

  IF(BUILD_LIBRA)
    ADD_DEFINITIONS(-DBUILD_LIBRA)
  ENDIF()

  IF(BUILD_CONFETTI)
    ADD_DEFINITIONS(-DBUILD_CONFETTI)
  ENDIF()

  IF(BUILD_FETALBRAIN)
    ADD_DEFINITIONS(-DBUILD_FETALBRAIN)
  ENDIF()

  IF(BUILD_GEODESICTRAINING)
    ADD_DEFINITIONS(-DBUILD_GEODESICTRAINING)
  ENDIF()

  # set stand-alone applications and libs directories
  SET_PROPERTY( GLOBAL PROPERTY USE_FOLDERS ON )
  SET( CaPTk_CLIs_Folder "CaPTk_CLIApps")
  SET( CaPTk_Libs_Folder "CaPTk_Libs")

  ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src)
  SET(INDIVIDUAL_APPLICATIONS_DIR ${PROJECT_SOURCE_DIR}/src/applications/individualApps)
  
  # include dependencies
  # add source code include directories
  INCLUDE_DIRECTORIES(
    ${CACHED_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/src/
    ${PROJECT_SOURCE_DIR}/src/applications/
    ${PROJECT_SOURCE_DIR}/src/applications/common_includes/
    ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/src/depends/
    ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/src/
    ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction_includes/
    ${PROJECT_SOURCE_DIR}/src/applications/GeodesicTraining/
    ${PROJECT_SOURCE_DIR}/src/applications/GeodesicTraining/include/
    ${PROJECT_SOURCE_DIR}/src/cbica_toolkit/src/
    ${PROJECT_SOURCE_DIR}/src/cbica_toolkit/src/IO/inc
    ${PROJECT_SOURCE_DIR}/src/common/inc/
    ${PROJECT_SOURCE_DIR}/src/view/gui/
    ${GEODESIC_TRAINING_INCLUDE_DIRS}
    ${APPLICATION_INCLUDES}
    ${CACHED_INCLUDE_DIRS}
    ${OpenMP_INCLUDE_DIRS}
  )
  
  IF( BUILD_CONFETTI )
    IF (WIN32)
      # do nothing for windows since the pre-compiled binaries will be used
    ELSE()

      SET(ConfettiSrcs ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/main.cpp
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confetti.cpp)
      SET(ConfettiHeaders
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confetti.h
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/dtiUtils.h
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/multinomialModel.h)

      #CAPTK_ADD_EXECUTABLE(${STANDALONE_APP_NAME} "${ConfettiSrcs} ${ConfettiHeaders}" "${OpenCV_LIBS} ${ITK_LIBRARIES}")
      
      ADD_EXECUTABLE(Confetti ${ConfettiSrcs} ${ConfettiHeaders})
      TARGET_LINK_LIBRARIES(Confetti ${OpenCV_LIBS} ${ITK_LIBRARIES})
      SET_TARGET_PROPERTIES( Confetti PROPERTIES FOLDER "${CaPTk_CLIApps}" )

      ADD_LIBRARY(ConfettiCore SHARED
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confetti.cpp
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confetti.h
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/multinomialModel.h
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/dtiUtils.h
        ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confettiApi.cpp
      )
      TARGET_LINK_LIBRARIES(ConfettiCore ${OpenCV_LIBS} ${ITK_LIBRARIES})

      # install all targets and scripts

      IF (APPLE) 
        INSTALL( TARGETS Confetti
        BUNDLE DESTINATION .
        RUNTIME DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
        LIBRARY DESTINATION ${EXE_NAME}.app/Contents/Resources/lib
        CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        )

        INSTALL( TARGETS ConfettiCore
        DESTINATION ${EXE_NAME}.app/Contents/Resources/lib
        )

        INSTALL( FILES ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiCore.py
          DESTINATION ${EXE_NAME}.app/Contents/Resources/libexec
        )

        IF( NOT CAPTK_CLI_MODE )
          INSTALL( FILES ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiGUI.py
            DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
            #CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
          )
        ENDIF()
      ELSE()
        INSTALL( TARGETS Confetti
        BUNDLE DESTINATION .
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        )

        INSTALL( TARGETS ConfettiCore
        DESTINATION lib
        )
        INSTALL( FILES ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiCore.py
          DESTINATION libexec
        )

        IF( NOT CAPTK_CLI_MODE )
          INSTALL( FILES ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiGUI.py
            DESTINATION bin
            #CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
          )
        ENDIF()
      ENDIF()

    ENDIF()
  ENDIF()

  ADD_DEFINITIONS(-DDATA_DIR="${DATA_DIR}")
  ADD_DEFINITIONS(-DICON_DIR="${DATA_DIR}/icons")
  ADD_DEFINITIONS(-DCAPTK_APP_LIST="${CAPTK_APP_LIST}")
  ADD_DEFINITIONS(-DPREPROCESS_ALGOS="${PREPROCESSING_ALGOS}")
  ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_GUI="${CAPTK_APP_LIST_PY_GUI}")
  ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_CLI="${CAPTK_APP_LIST_PY_CLI}")

  # keep source files structured as folders for clarity
  SOURCE_GROUP( GUI FILES ${GUI_SOURCES_UPDATED} )
  SOURCE_GROUP( CBICA-TK FILES ${CBICA-TK_SOURCES} )
  SOURCE_GROUP( APPLICATIONS FILES ${APPLICATION_SOURCES} )
  SOURCE_GROUP( APPLICATIONS_Dependencies FILES ${APPLICATION_DEPENDS} )

  #SET( LIBNAME_GUI "${PROJECT_NAME}_GUI" )
  SET( LIBNAME_Applications "${PROJECT_NAME}_Applications" )
  SET( LIBNAME_FeatureExtractor "${PROJECT_NAME}_FeatureExtractor" )
  SET( LIBNAME_CBICATK "${PROJECT_NAME}_CBICATK" )

  # add all dependenciesO of ${application} to compilation
  FILE( GLOB 
    FEATURE_EXTRACTION_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction_includes/*.h"
    "${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction_includes/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction_includes/*.hxx"
  )

  ADD_LIBRARY(
    ${LIBNAME_FeatureExtractor}
    ${FEATURE_EXTRACTION_DEPENDS}
    ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/src/FeatureExtraction.h
    ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/src/FeatureExtraction.hxx
    ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/src/depends/LBPFeatures2D.h
    ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/src/depends/LBPFeatures2D.cpp
    #${CBICA-TK_SOURCES}
  )

  TARGET_LINK_LIBRARIES(
    ${LIBNAME_FeatureExtractor}
    #${DEPENDENT_LIBS}
    ${LIBNAME_CBICATK}
    ${ITK_LIBRARIES}
    ${OpenCV_LIBRARIES}
  )

  ADD_LIBRARY(
    ${LIBNAME_Applications}
    ${APPLICATION_SOURCES}
    ${APPLICATION_DEPENDS}
    #${CBICA-TK_SOURCES}
  )

  TARGET_LINK_LIBRARIES(
    ${LIBNAME_Applications}
    ${DEPENDENT_LIBS}
    ${LIBNAME_CBICATK}
    ${ITK_LIBRARIES}
    ${OpenCV_LIBRARIES}
  )

  ADD_DEPENDENCIES( ${LIBNAME_Applications} ${LIBNAME_GUI} ${LIBNAME_FeatureExtractor} )

  SET_TARGET_PROPERTIES( ${LIBNAME_CBICATK} PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
  SET_TARGET_PROPERTIES( ${LIBNAME_Applications} PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
  SET_TARGET_PROPERTIES( ${LIBNAME_FeatureExtractor} PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
  SET_TARGET_PROPERTIES( CaPTk_GUI PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
  SET_TARGET_PROPERTIES( CaPTk_Common PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
  SET_TARGET_PROPERTIES( GeodesicTrainingCaPTkLib PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
  SET_TARGET_PROPERTIES( GeodesicTrainingComputeLib PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )

  
  SET( CAPTK_LIBRARIES 
    ${LIBNAME_CBICATK} 
    ${LIBNAME_Applications} 
    ${LIBNAME_FeatureExtractor} 
  )

  # find stand-alone CMD applications
  FILE( GLOB STANDALONE_APPS
    "${PROJECT_SOURCE_DIR}/src/applications/*.cxx"
  )

  #MESSAGE( STATUS "[debug]STANDALONE_APPS == ${STANDALONE_APPS}")

  #SET(STANDALONE_APPS_LIST "")

  FOREACH(application ${STANDALONE_APPS})

    GET_FILENAME_COMPONENT( STANDALONE_APP_NAME ${application} NAME_WE )
    
    IF( NOT (${STANDALONE_APP_NAME} STREQUAL "ImagingSubtypePredictor" ) )
    
        CAPTK_ADD_EXECUTABLE(${STANDALONE_APP_NAME} ${PROJECT_SOURCE_DIR}/src/applications/${STANDALONE_APP_NAME}.cxx "${CAPTK_LIBRARIES}" )
    
    ENDIF()
  ENDFOREACH()

  IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install" CACHE PATH "Default Install Path" FORCE )
  ENDIF()

  ## Build test(s) only if CAPTK_PACKAGE_PROJECT is OFF
  #IF( NOT CAPTK_PACKAGE_PROJECT )
    SET( BUILD_TESTING TRUE CACHE BOOL "Build testing" )
    IF( BUILD_TESTING )
      ENABLE_TESTING()
      ADD_SUBDIRECTORY( testing )
      #INCLUDE(CTest)
    ENDIF( BUILD_TESTING )
  #ENDIF()

  # installation
  IF (UNIX)
    IF(NOT LIBEXEC_INSTALL_DIR)
      SET(LIBEXEC_INSTALL_DIR "./libexec")
    ENDIF()
  ENDIF()

  ## symlink creation for UNIX machines - currently only works for DEB or RPM
  # IF (UNIX OR APPLE)
  #   INSTALL(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/bin/${EXE_NAME} ${CMAKE_INSTALL_PREFIX}/bin/captk)")
  # ENDIF()

  SET(META_FILES_TO_INCLUDE "") # used for Apple package
  #IF( BUILD_DOCUMENTATION )
    FILE(GLOB ALL_DOC_FILES_BASE "${PROJECT_BINARY_DIR}/docs_sources/html/*.*")
    FOREACH(docFiles ${ALL_DOC_FILES_BASE})
      IF (APPLE)
        SET_SOURCE_FILES_PROPERTIES(${docFiles} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/share/doc)
        SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${docFiles})
      ELSE()
        INSTALL(FILES "${docFiles}" DESTINATION share/doc)
      ENDIF()
    ENDFOREACH()
    #INSTALL( DIRECTORY ${PROJECT_BINARY_DIR}/docs_sources/html DESTINATION share/doc )
  #ENDIF()

  # ensure all licenses are in the package
  FILE( GLOB_RECURSE ALL_LICENSES "${PROJECT_SOURCE_DIR}/licenses/*.txt" )
  IF (APPLE)
    SET_SOURCE_FILES_PROPERTIES(${ALL_LICENSES} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/license)
    SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${ALL_LICENSES})

  ELSE()
    FOREACH( license ${ALL_LICENSES} )
      INSTALL( FILES "${license}" DESTINATION licenses )
    ENDFOREACH()
  ENDIF()

  #INSTALL(FILES "${README_FILE}" DESTINATION ".") # not needed any more

  # ensure all icons are in the package
  FILE( GLOB_RECURSE ALL_ICONS "${DATA_DIR}/icons/*.png" )
  IF (APPLE)
    SET_SOURCE_FILES_PROPERTIES(${ALL_ICONS} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/icons)
    SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${ALL_ICONS})

  ELSE()
    FOREACH( icon ${ALL_ICONS} )
      INSTALL(FILES "${icon}" DESTINATION data/icons )
    ENDFOREACH()
  ENDIF()

  SET( CAPTK_STYLESHEET_FILE "captk.qss" )
  ADD_DEFINITIONS( -DCAPTK_STYLESHEET_FILE="${CAPTK_STYLESHEET_FILE}" )
  #INSTALL(FILES "${DATA_DIR}/FetalModel.xml" DESTINATION data)
  IF (WIN32)
    # Can't be used with Linuxdeployqt or macdeployqt
    INSTALL(FILES "${DATA_DIR}/${CAPTK_STYLESHEET_FILE}" DESTINATION etc)
    INSTALL(FILES "${DATA_DIR}/qt.conf" DESTINATION bin)
  ELSEIF (APPLE)
    SET_SOURCE_FILES_PROPERTIES("${DATA_DIR}/${CAPTK_STYLESHEET_FILE}" PROPERTIES MACOSX_PACKAGE_LOCATION Resources/etc)
    SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} "${DATA_DIR}/${CAPTK_STYLESHEET_FILE}")

  ELSE ()
    INSTALL(FILES "${DATA_DIR}/${CAPTK_STYLESHEET_FILE}" DESTINATION etc)
    # Linux needs the applications folder to be able to create a binary with linuxdeployqt
    INSTALL(DIRECTORY "${DATA_DIR}/applications" DESTINATION share)
  ENDIF()

  IF (APPLE)
    SET( FLAIR_FILE "${DATA_DIR}/AAAC0_flair_pp_shrunk.nii.gz")
    SET_SOURCE_FILES_PROPERTIES( ${FLAIR_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data)

    SET( FLAIR_TESTTUMOR_FILE "${DATA_DIR}/AAAC0_flair_pp_shrunk_testTumor.nii.gz")
    SET_SOURCE_FILES_PROPERTIES( ${FLAIR_TESTTUMOR_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data)

    FILE(GLOB survival_FILES "${DATA_DIR}/survival/*.*")
    SET_SOURCE_FILES_PROPERTIES(${survival_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/survival)
    
    FILE(GLOB egfrv3_FILES "${DATA_DIR}/egfrv3/*.*")
    SET_SOURCE_FILES_PROPERTIES(${egfrv3_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/egfrv3)
    
    FILE(GLOB egfrv3_FILES "${DATA_DIR}/molecular/*.*")
    SET_SOURCE_FILES_PROPERTIES(${egfrv3_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/molecular)
    
    FILE(GLOB sri24_FILES "${DATA_DIR}/sri24/*.*")
    SET_SOURCE_FILES_PROPERTIES(${sri24_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/sri24)

    FILE(GLOB DOWNLOAD_LINKS "${DATA_DIR}/links.yaml")
    SET_SOURCE_FILES_PROPERTIES(${DOWNLOAD_LINKS} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data)

    SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${FLAIR_FILE} ${FLAIR_TESTTUMOR_FILE} ${survival_FILES} ${egfrv3_FILES} ${sri24_FILES} ${DOWNLOAD_LINKS})

  ELSE ()
    INSTALL(FILES "${DATA_DIR}/AAAC0_flair_pp_shrunk.nii.gz" DESTINATION data)
    INSTALL(FILES "${DATA_DIR}/AAAC0_flair_pp_shrunk_testTumor.nii.gz" DESTINATION data)
    INSTALL(FILES "${DATA_DIR}/links.yaml" DESTINATION data)
    INSTALL(DIRECTORY "${DATA_DIR}/survival" DESTINATION data PATTERN ".svn" EXCLUDE)
    INSTALL(DIRECTORY "${DATA_DIR}/egfrv3" DESTINATION data PATTERN ".svn" EXCLUDE)
    INSTALL(DIRECTORY "${DATA_DIR}/molecular" DESTINATION data PATTERN ".svn" EXCLUDE)
    INSTALL(DIRECTORY "${DATA_DIR}/sri24" DESTINATION data PATTERN ".svn" EXCLUDE)
  ENDIF()

  # install DeepMedic files
  IF (APPLE)
    FILE(GLOB_RECURSE DEEPMEDIC_CONFIGFILES "${DATA_DIR}/deepMedic/configFiles/*")
    FILE(GLOB_RECURSE DEEPMEDIC_SAVEDMODELS_BRAIN_TUMOR "${DATA_DIR}/deepMedic/saved_models/brainTumorSegmentation/*")
    FILE(GLOB_RECURSE DEEPMEDIC_SAVEDMODELS_SKULL "${DATA_DIR}/deepMedic/saved_models/skullStripping/*")
    
    SET_SOURCE_FILES_PROPERTIES( ${DEEPMEDIC_CONFIGFILES} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/deepMedic/configFiles )
    SET_SOURCE_FILES_PROPERTIES( ${DEEPMEDIC_SAVEDMODELS_BRAIN_TUMOR} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/deepMedic/saved_models/brainTumorSegmentation )
    SET_SOURCE_FILES_PROPERTIES( ${DEEPMEDIC_SAVEDMODELS_SKULL} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/deepMedic/saved_models/skullStripping )
    
    SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${DEEPMEDIC_CONFIGFILES})
    SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${DEEPMEDIC_SAVEDMODELS_BRAIN_TUMOR})
    SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${DEEPMEDIC_SAVEDMODELS_SKULL})
    
  ELSE()
    INSTALL(DIRECTORY "${DATA_DIR}/deepMedic/configFiles" DESTINATION data/deepMedic)
    INSTALL(DIRECTORY "${DATA_DIR}/deepMedic/saved_models" DESTINATION data/deepMedic)
  ENDIF()
  
  # install CWL files if they aren't getting generated
  #IF(CAPTK_DISABLE_CWL)
	
    FILE(GLOB_RECURSE CWL_FILES "${DATA_DIR}/cwlFiles/*")
    FOREACH(cwlFile ${CWL_FILES})
      IF (APPLE)
        SET_SOURCE_FILES_PROPERTIES(${cwlFile} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/etc/cwlDefinitions)
        SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${cwlFile})

        INSTALL(FILES "${cwlFile}" DESTINATION ${EXE_NAME}.app/Contents/Resources/etc/cwlDefinitions)

      ELSE()
        INSTALL(FILES "${cwlFile}" DESTINATION etc/cwlDefinitions)
      ENDIF()
    ENDFOREACH()
  
  #ENDIF()

  IF(WIN32)
    FILE(GLOB_RECURSE ALL_DEPENDENT_LIBS "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/*.dll")
  ELSE(WIN32)
    FILE(GLOB_RECURSE ALL_DEPENDENT_LIBS "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/*.so")
  ENDIF()

  FOREACH(dependent_lib ${ALL_DEPENDENT_LIBS})
    IF (APPLE)
      SET_SOURCE_FILES_PROPERTIES(${dependent_lib} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/bin)
      SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${dependent_lib})

      INSTALL(FILES "${dependent_lib}" DESTINATION ${EXE_NAME}.app/Contents/Resources/bin)

    ELSE()
      INSTALL(FILES "${dependent_lib}" DESTINATION bin)
    ENDIF()
  ENDFOREACH()

  #IF( WIN32 )
  #  INSTALL(DIRECTORY "${DATA_DIR}/GnuWin32" DESTINATION data PATTERN ".svn" EXCLUDE)
  #ENDIF()

  FILE(GLOB ALL_DOC_FILES_BASE "${PROJECT_BINARY_DIR}/docs_sources/html/*.*")
  FOREACH(docFiles ${ALL_DOC_FILES_BASE})
    IF (APPLE)
      SET_SOURCE_FILES_PROPERTIES(${docFiles} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/share/doc)
      SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${docFiles})

    ELSE()
      INSTALL(FILES "${docFiles}" DESTINATION share/doc)
    ENDIF()
  ENDFOREACH()

  set(OS_BUNDLE MACOSX_BUNDLE)

  IF( NOT CAPTK_CLI_MODE )

    # define the executable
    IF(WIN32)
      # ensures no console pops up when launching the program
      ADD_EXECUTABLE(${EXE_NAME} WIN32
        ${EXE_SOURCES}
        #${YAMLCPP_Headers}
        #${YAMLCPP_Sources}
        ${DATA_DIR}/icons/application/windows/captk.rc
      )
    ELSEIF(APPLE)
      SET (ICON_FILE "${DATA_DIR}/icons/application/apple/CaPTk.icns")
      SET_SOURCE_FILES_PROPERTIES(${ICON_FILE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

      FOREACH(CAPTK_MACOSX_BUNDLE_FILE ${CAPTK_MACOSX_BUNDLE_FILES_INSTALL_TARGETS}) 
          # MESSAGE(STATUS ${CAPTK_MACOSX_BUNDLE_FILES_INSTALL_TARGETS}) 
          SET_SOURCE_FILES_PROPERTIES(${CAPTK_MACOSX_BUNDLE_FILE} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE PROPERTIES MACOSX_PACKAGE_LOCATION Resources/bin )

          SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${CAPTK_MACOSX_BUNDLE_FILE})
      ENDFOREACH()

      FOREACH(FT_META_FILES_TO_INCLUDE ${FEATURE_EXTRACTION_META_FILES_TO_INCLUDE})
        IF (${FT_META_FILES_TO_INCLUDE} STREQUAL "${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/data/batchMode/batch_featureExtraction.csv")
          SET_SOURCE_FILES_PROPERTIES(${FT_META_FILES_TO_INCLUDE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/share/featureExtractionBatch )

          SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${FT_META_FILES_TO_INCLUDE})
        ELSE() 
          SET_SOURCE_FILES_PROPERTIES(${FT_META_FILES_TO_INCLUDE} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/data/features )

          SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${FT_META_FILES_TO_INCLUDE}) 
        ENDIF()        
      ENDFOREACH()

      IF (NOT ICON_NAME)
        SET(ICON_NAME CaPTk)
      ENDIF()

      #   SET_SOURCE_FILES_PROPERTIES(Confetti PROPERTIES MACOSX_PACKAGE_LOCATION Resources/bin)
      #   SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} Confetti)

      # # SET_SOURCE_FILES_PROPERTIES(${ConfettiCore} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/lib)
      # # SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${ConfettiCore})

      SET_SOURCE_FILES_PROPERTIES("${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiGUI.py" PROPERTIES MACOSX_PACKAGE_LOCATION Resources/bin)
      SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} "${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiGUI.py")

      SET_SOURCE_FILES_PROPERTIES("${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiCore.py" PROPERTIES MACOSX_PACKAGE_LOCATION Resources/libexec)
      SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} "${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/pyGUI/ConfettiCore.py")
      
      # FOREACH(application ${STANDALONE_APPS})

      #   GET_FILENAME_COMPONENT( STANDALONE_APP_NAME ${application} NAME_WE )
    
      #   IF( NOT (${STANDALONE_APP_NAME} STREQUAL "ImagingSubtypePredictor" ) )
      #       SET_SOURCE_FILES_PROPERTIES(${STANDALONE_APP_NAME} PROPERTIES MACOSX_PACKAGE_LOCATION Resources/bin)
      #       SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} ${STANDALONE_APP_NAME})
      #   ENDIF()
      # ENDFOREACH()

      #   SET_SOURCE_FILES_PROPERTIES(CaPTk_Utilities PROPERTIES MACOSX_PACKAGE_LOCATION Resources/bin)
      #   SET(META_FILES_TO_INCLUDE ${META_FILES_TO_INCLUDE} CaPTk_Utilities)

      # Identify MacOS bundle
      SET(MACOSX_BUNDLE_BUNDLE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}")
      SET(MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION})
      SET(MACOSX_BUNDLE_LONG_VERSION_STRING ${PROJECT_VERSION})
      SET(MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION})
      SET(MACOSX_BUNDLE_COPYRIGHT ${COPYRIGHT})
      SET(MACOSX_BUNDLE_GUI_IDENTIFIER ${IDENTIFIER})
      SET(MACOSX_BUNDLE_ICON_FILE ${ICON_NAME})

      ADD_EXECUTABLE(${EXE_NAME}
        ${OS_BUNDLE}
        ${EXE_SOURCES}
        ${ICON_FILE}
        ${META_FILES_TO_INCLUDE}
      )
    ELSE()
      ADD_EXECUTABLE(${EXE_NAME}
        ${EXE_SOURCES}
        #${YAMLCPP_Headers}
        #${YAMLCPP_Sources}
      )
    ENDIF()

    TARGET_LINK_LIBRARIES( ${EXE_NAME}
      ${DEPENDENT_LIBS}
      CaPTk_GUI
      CaPTk_Common
      GeodesicTrainingComputeLib
      GeodesicTrainingCaPTkLib
      ${LIBNAME_Applications}
      ${LIBNAME_FeatureExtractor}
      ${LIBNAME_CBICATK}
    )

    ADD_DEPENDENCIES( ${EXE_NAME}
      CaPTk_GUI
      CaPTk_Common
      ${LIBNAME_Applications}
      ${LIBNAME_FeatureExtractor}
      ${LIBNAME_CBICATK}
    )

    # set properties for the executable
    # SET_TARGET_PROPERTIES( ${EXE_NAME}
    #   PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE
    #   OUTPUT_NAME ${EXE_NAME}
    # )

  ENDIF()

  IF (APPLE)
    INSTALL( TARGETS ${EXE_NAME}
      BUNDLE DESTINATION .
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

    # INSTALL( FILES "${PROJECT_BINARY_DIR}/${EXE_NAME}.app/Contents/MacOS/${EXE_NAME}"
    #   DESTINATION ${EXE_NAME}.app/Contents/MacOS
    #   PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    # )
    
  ELSE()
    INSTALL( TARGETS ${EXE_NAME}
      BUNDLE DESTINATION .
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

  ENDIF()

  IF(CAPTK_PACKAGE_PROJECT)
  
    SET(DEPENDENCY_DIRS ${VTK_LIBRARY_DIRS}/${CMAKE_CONFIGURATION_TYPES}/ ${QT_BINARY_DIR})
    IF(APPLE)
      # SET(EXECUTABLE "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/${EXE_NAME}.app")
    ELSEIF(WIN32)
      SET(EXECUTABLE "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/${EXE_NAME}.exe")
    ELSE()
      SET(EXECUTABLE "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/${EXE_NAME}")
    ENDIF()

    SET(CPACK_PACKAGE_EXECUTABLES "${EXECUTABLE}" "${EXE_NAME}")

    IF(WIN32)
      SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}_${${PROJECT_NAME}_Variant}\\\\${PROJECT_VERSION}")
    ELSE()
      SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}_${${PROJECT_NAME}_Variant}/${PROJECT_VERSION}")
    ENDIF()

    IF( WIN32 )
      INSTALL(FILES "${DATA_DIR}/icons/application/windows/captk.ico" DESTINATION share/icons)
      INSTALL(FILES "${DATA_DIR}/icons/application/windows/captk.rc" DESTINATION data/icons)
    ELSEIF (APPLE)
      
    ELSE()
      INSTALL( DIRECTORY "${DATA_DIR}/icons/application/linux/hicolor" DESTINATION share/icons)
      # Not in data dir on linux for some reason
      # INSTALL(FILES "${DATA_DIR}/icons/captk.rc" DESTINATION data/icons)
    ENDIF()

    # CPack settings
    IF(WIN32)
      SET(CPACK_GENERATOR NSIS64)
      SET(CPACK_SOURCE_GENERATOR "ZIP")
      SET(CPACK_NSIS_DISPLAY_NAME "${PROJECT_NAME}")
      SET(CPACK_NSIS_PACKAGE_NAME "${PROJECT_NAME}")
      SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}_${${PROJECT_NAME}_Variant}\\\\${PROJECT_VERSION}")
      SET(CPACK_NSIS_INSTALL_ROOT "C:")
      SET(CPACK_NSIS_DEFINES
        "RequestExecutionLevel user
        ManifestDPIAware true"
      )
      #SET(CPACK_NSIS_DEFINES "RequestExecutionLevel user ManifestDPIAware true")
      SET(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
      SET(CPACK_NSIS_MUI_ICON "${DATA_DIR}/icons/application/windows/captk.ico")
      SET(CPACK_NSIS_MUI_FINISHPAGE_RUN "${EXE_NAME}.exe")
      SET(CPACK_NSIS_HELP_LINK "https:\\\\\\\\www.med.upenn.edu\\\\cbica\\\\captk\\\\")
      SET(CPACK_NSIS_CONTACT "${PROJECT_CONTACT}")
      SET(CPACK_NSIS_URL_INFO_ABOUT "https:\\\\\\\\www.med.upenn.edu\\\\cbica\\\\captk\\\\")
      #SET(CPACK_NSIS_MODIFY_PATH "ON") # desktop shortcut doesn't work, for some reason
      # create extra shortcuts and respective uninstall targets
      SET(CPACK_NSIS_CREATE_ICONS_EXTRA "CreateShortCut '\$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\${EXE_NAME}.lnk' '\$INSTDIR\\\\bin\\\\${EXE_NAME}.exe'")
      #SET(CPACK_NSIS_CREATE_ICONS "CreateShortCut '$DESKTOP\\\\${EXE_NAME}.lnk' '$INSTDIR\\\\bin\\\\${EXE_NAME}.exe'")
      #SET(CPACK_NSIS_CREATE_ICONS_EXTRA "  CreateShortCut '$INSTDIR\\\\${EXE_NAME}.lnk' '$INSTDIR\\\\bin\\\\${EXE_NAME}.exe'")
      SET(CPACK_NSIS_CREATE_ICONS "CreateShortCut '\$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\README.lnk' '\$INSTDIR\\\\share\\\\doc\\\\index.html'")
      #SET(CPACK_NSIS_DELETE_ICONS_EXTRA "  Delete '$SMPROGRAMS\\\\$START_MENU\\\\${EXE_NAME}.lnk'")
      SET(CPACK_NSIS_DELETE_ICONS 
        "Delete '$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\${EXE_NAME}.lnk'
        Delete '$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\README.lnk'"
      )
      #SET(CPACK_NSIS_DELETE_ICONS "Delete '$DESKTOP\\\\${EXE_NAME}.lnk'")
      #SET(CPACK_NSIS_DELETE_ICONS_EXTRA "Delete '$INSTDIR\\\\${EXE_NAME}.lnk'")
      #SET(CPACK_NSIS_DELETE_ICONS "Delete '$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\README.lnk'")
      SET(CPACK_NSIS_INSTALLED_ICON_NAME "${EXE_NAME}")
      SET(CPACK_SYSTEM_NAME "Win64")
    ELSEIF(APPLE)
      SET(CPACK_GENERATOR DragNDrop )
      SET(CPACK_DMG_FORMAT "UDBZ")
      SET(CPACK_DMG_VOLUME_NAME "${PROJECT_NAME}_${PROJECT_VERSION}")
      SET(CPACK_BUNDLE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}")
      SET(CPACK_SYSTEM_NAME "OSX")
      # SET(CPACK_STRIP_FILES "${EXE_NAME}")
      # SET(CPACK_STRIP_FILES "${EXE_NAME}.app/Contents/MacOS/${EXE_NAME}")
    ENDIF()
    # common
    SET(CPACK_SOURCE_STRIP_FILES TRUE )
    SET(CPACK_PACKAGE_EXECUTABLES "${EXE_NAME}" "${EXE_NAME}")
    SET(CPACK_CREATE_DESKTOP_LINKS "${EXE_NAME}" )
    SET(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}")
    IF(WIN32)
      SET(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}_Installer")
    ENDIF()
    SET(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
    SET(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
    SET(CPACK_PACKAGE_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
    SET(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
    SET(CPACK_PACKAGE_VENDOR "${PROJECT_VENDOR}")
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CaPTk is a collection of algorithms from CBICA")
    SET(CPACK_SOURCE_IGNORE_FILES "\\\\.svn/;^${PROJECT_SOURCE_DIR}/doc/")
    SET(CPACK_SOURCE_PACKAGE_FILE_NAME "captk-${SW_VER}-source" CACHE INTERNAL "tarball basename")
    SET(CPACK_RESOURCE_FILE_README "${README_FILE}")
    SET(CPACK_RESOURCE_FILE_LICENSE "${LICENSE_FILE}")
    SET(CPACK_COMPONENT_LIBRARIES_INSTALL_TYPES Developer Full)
    SET(CPACK_COMPONENT_HEADERS_INSTALL_TYPES Developer Full)
    SET(CPACK_COMPONENT_APPLICATIONS_INSTALL_TYPES Full)
    SET(CPACK_ALL_INSTALL_TYPES Full Developer)
    SET(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/bin/;{CMAKE_SOURCE_DIR}/build/;${CMAKE_SOURCE_DIR}/.svn/")

    SET(CMAKE_INSTALL_UCRT_LIBRARIES TRUE)
    SET(CMAKE_INSTALL_OPENMP_LIBRARIES TRUE)
    INCLUDE(InstallRequiredSystemLibraries)

    IF(WIN32)
      INCLUDE(Windeployqt)
      windeployqt(${EXE_NAME} bin)
    ELSEIF(APPLE)
      INCLUDE(ChmodFile)
      chmodFile(${EXE_NAME} dcm2nii)
      
      INCLUDE(Macdeployqt)
      macdeployqt(${EXE_NAME})
      
      INCLUDE(MvFile)
      mvFile(${EXE_NAME} itksnap ITK-SNAP)
      chmodFile(${EXE_NAME} "ITK-SNAP.app/Contents/MacOS/ITK-SNAP")

      INCLUDE(ChangeDylib)

      FILE( GLOB EXECS
        "${DATA_DIR}/cwlFiles/*"
      )

      FOREACH(EXEC ${EXECS})

        GET_FILENAME_COMPONENT( EXEC_IN_BIN ${EXEC} NAME_WE )

        # GET_TARGET_PROPERTY(TARGET_TYPE ${EXEC_IN_BIN} TYPE)
        IF (NOT EXEC_IN_BIN STREQUAL "BreastTexturePipeline")
          # MESSAGE( STATUS "[debug]EXEC_IN_BIN == ${EXEC_IN_BIN}")
          changeLibOmpDylib(${EXE_NAME} ${EXEC_IN_BIN})
        ENDIF ()
        
      ENDFOREACH()

    ENDIF()

    INCLUDE(CPack)

    #IF(EXISTS ${EXE_NAME}) # don't do anything if the EXE is not found - doing this instead of "install_qt4_exe" since it works
    #  # this will no longer be needed with qt5 since it comes with its own packaging tool
    #  #IF(APPLE OR WIN32)
    #    #INCLUDE(DeployQt4)
    #    #INCLUDE(BundleUtilities)
    #    #MESSAGE(STATUS "Qt libs: ${QT_QTCORE_LIBRARY};${QT_QTGUI_LIBRARY}")
    #    #MESSAGE(STATUS "Qt lib dir: ${QT_LIBRARY_DIR}")
    #    #MESSAGE(STATUS "Qt bin dir: ${QT_BINARY_DIR}")
    #    #MESSAGE(STATUS "VTK bin dir: ${VTK_LIBRARY_DIRS}/${CMAKE_CONFIGURATION_TYPES}/")
    #    #FIXUP_BUNDLE("${EXECUTABLE}" "" "${DEPENDENCY_DIRS}")
    #    #FIXUP_QT4_EXECUTABLE("${EXECUTABLE}" "" "" "${QT_LIBRARY_DIR}" "" "")
    #    #SET(qtconf_dest_dir /Contents/Resources)
    #    #INSTALL(CODE
    #    #  "
    #    #  FILE(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"\")
    #    #  "
    #    #  COMPONENT Runtime
    #    #)
    #    IF(WIN32)
    #      INSTALL( CODE
    #        "
    #        INCLUDE(BundleUtilities)
    #        FIXUP_BUNDLE( \"${EXE_NAME}\" \"\" \"${QT_BINARY_DIR}\" )
    #        "
    #        Component Runtime
    #      )
    #    #INSTALL_QT4_EXECUTABLE("${EXECUTABLE}" "" "${QT_LIBS}" "${QT_LIBRARY_DIR}" "" "" "")
    #    ENDIF()
    #  #ENDIF()
    #ENDIF()

  ENDIF(CAPTK_PACKAGE_PROJECT)

  # Execute a post install step. Make sure this is the last install command ran!
  IF( CMAKE_HOST_UNIX )
    INSTALL(CODE "EXECUTE_PROCESS(COMMAND ln -sf CaPTk captk WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin)")
    # INSTALL(CODE "EXECUTE_PROCESS(COMMAND chmod 755 matlab_helper WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/libexec/MCR/v84/bin/glnxa64/)")

    # IF(EXISTS ${PROJECT_SOURCE_DIR}/src/applications/individualApps/libra/MCR/v84/bin/glnxa64/matlab_helper)
    #   INSTALL(CODE "MESSAGE(\"Trying to install matlab_helper with correct permissions.\")")
    #   INSTALL(FILES ${PROJECT_SOURCE_DIR}/src/applications/individualApps/libra/MCR/v84/bin/glnxa64/matlab_helper
    #     DESTINATION libexec/MCR/v84/bin/glnxa64/
    #     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    #   )
    # ELSEIF(EXISTS ${PROJECT_BINARY_DIR}/externalApps/binaries_linux/libra/MCR/v84/bin/glnxa64/matlab_helper)
    #   INSTALL(CODE "MESSAGE(\"Trying to install matlab_helper with correct permissions.\")")
    #   INSTALL(FILES ${PROJECT_BINARY_DIR}/externalApps/binaries_linux/libra/MCR/v84/bin/glnxa64/matlab_helper
    #     DESTINATION libexec/MCR/v84/bin/glnxa64/
    #     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    #   )
    # ENDIF()

  ENDIF()

ENDIF( NOT BUILD_DOCUMENTATION_ONLY )
