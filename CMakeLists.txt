# parent CMakeLists file

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)
IF(COMMAND cmake_policy)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0020 NEW)
  CMAKE_POLICY(SET CMP0010 NEW)
  #CMAKE_POLICY(SET CMP0048 NEW) # only with CMake 3
ENDIF(COMMAND cmake_policy)

# set project name. The exe name is taken from this
PROJECT( CaPTk )
SET( EXE_NAME "${PROJECT_NAME}" )
SET( ${PROJECT_NAME}_Variant "Full" ) # the particular variant of CaPTk (Full/Neuro/Torso/Viewer)
SET (PROJECT_VERSION_MAJOR 1 )
SET (PROJECT_VERSION_MINOR 5 )
SET (PROJECT_VERSION_PATCH 0 )
SET( PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
IF(${PROJECT_NAME}_Variant MATCHES "Full" )
  SET( PROJECT_NAME_EXTENDED "Cancer Imaging Phenomics Toolkit (${PROJECT_NAME})")
ELSE()  
  SET( PROJECT_NAME_EXTENDED "Cancer Imaging Phenomics Toolkit (${PROJECT_NAME}_${${PROJECT_NAME}_Variant})" )
ENDIF()
SET( PROJECT_VENDOR "CBICA - UPenn" )
SET( LICENSE_FILE "${PROJECT_SOURCE_DIR}/licenses/Combined.txt" )
SET( README_FILE  "${PROJECT_SOURCE_DIR}/README.txt" )
SET( PROJECT_DOMAIN "http://www.med.upenn.edu/cbica/captk/" )
SET( PROJECT_CONTACT "software@cbica.upenn.edu" )
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

# Provide build configuration option
SET( CAPTK_BUILD_OPTION "PACKAGE_PROJECT" CACHE STRING "Toggle between packaging mode (optimized code) or developer mode (provides debugging info)" )

SET_PROPERTY( CACHE CAPTK_BUILD_OPTION PROPERTY STRINGS PACKAGE_PROJECT DEVELOPER_MODE )

IF( CAPTK_BUILD_OPTION STREQUAL PACKAGE_PROJECT)

  SET( CAPTK_PACKAGE_PROJECT ON )
  SET( CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Default configuration types" FORCE )
  ADD_DEFINITIONS(-DCAPTK_PACKAGE_PROJECT)

  # add a few compiler options for GCC and MSVC
  IF( CMAKE_COMPILER_IS_GNUCXX )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function" )
  ENDIF() 
  
  IF(MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
  ENDIF(MSVC)

ELSE()

  SET( CAPTK_PACKAGE_PROJECT OFF )
  SET( CMAKE_CONFIGURATION_TYPES "Debug" CACHE STRING "Default configuration types" FORCE )

ENDIF()

MARK_AS_ADVANCED(CAPTK_PACKAGE_PROJECT)
MARK_AS_ADVANCED(CAPTK_BUILD_OPTION)

SET( CMAKE_COLOR_MAKEFILE ON )

# add c++11 flag if it is GCC v<5
IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    INCLUDE( CheckCXXCompilerFlag )
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    IF( COMPILER_SUPPORTS_CXX11 )
    	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    elseif(COMPILER_SUPPORTS_CXX0X )
    	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else()
    	MESSAGE(ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
  ENDIF()
ENDIF() 

# option to compile only the viewer
SET( CAPTK_BUILD_CONSOLE_ONLY OFF CACHE BOOL "Do you want only the interactive capabilities of ${PROJECT_NAME}? Enabling this will remove all algorithms from compilation." )

IF (CAPTK_BUILD_CONSOLE_ONLY)
  ADD_DEFINITIONS(-DCAPTK_BUILD_CONSOLE_ONLY)
ENDIF()

MARK_AS_ADVANCED(CAPTK_BUILD_CONSOLE_ONLY)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)

SET( CAPTK_CLI_MODE OFF CACHE BOOL "Builds only CLI apps of CaPTk." )
ADD_DEFINITIONS(-DCAPTK_CLI_MODE)

# find dependencies START
FIND_PACKAGE(ITK REQUIRED)
INCLUDE("${ITK_USE_FILE}")

FIND_PACKAGE(OpenCV REQUIRED)
#MESSAGE( STATUS "OpenCV dirs:${OpenCV_INCLUDE_DIRS}")
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

FIND_PACKAGE(OpenMP REQUIRED)
SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

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

#INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/SuperBuild.cmake")

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(UseQt4)
INCLUDE(${QT_USE_FILE})
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTMAIN TRUE) 
SET(QT_USE_QTGUI TRUE)
ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)
LINK_DIRECTORIES(${QT_LIBRARY_DIR})

# link dependent libraries
SET( DEPENDENT_LIBS 
  Qt4::QtCore
  Qt4::QtGui
  #QVTK 
  ${ITK_LIBRARIES} 
  ${VTK_LIBRARIES}
  ${OpenCV_LIBRARIES}
)

# find dependencies END

# add source code include directories
INCLUDE_DIRECTORIES(
	${QT_INCLUDES}
	${QT_INCLUDE_DIR}
	${QT_QTGUI_INCLUDE_DIR}
	${QT_QTCORE_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/src/
  ${PROJECT_SOURCE_DIR}/src/applications/
  ${PROJECT_SOURCE_DIR}/src/applications/common_includes/
  ${PROJECT_SOURCE_DIR}/src/applications/WhiteStripe_includes/
  #${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction_includes/
  ${PROJECT_SOURCE_DIR}/src/cbica_toolkit/
  ${PROJECT_SOURCE_DIR}/src/gui/
  ${APPLICATION_INCLUDES}
)

SET( DATA_DIR "${PROJECT_SOURCE_DIR}/data/" )

### Building the documentation
## The following 4 lines are for generation of public documentation
SET( BUILD_DOCUMENTATION ON CACHE BOOL "Do you want to build Documentation (needs Doxygen)?" )
SET( BUILD_DOCUMENTATION_ONLY OFF CACHE BOOL "This doesn't build anything else in the project" )
MARK_AS_ADVANCED(BUILD_DOCUMENTATION)
MARK_AS_ADVANCED(BUILD_DOCUMENTATION_ONLY)
IF( BUILD_DOCUMENTATION_ONLY )
  SET( BUILD_DOCUMENTATION ON )
ENDIF()
IF( BUILD_DOCUMENTATION )
  ADD_SUBDIRECTORY( docs )
ENDIF ()

IF( NOT BUILD_DOCUMENTATION_ONLY )
# define build options for applications and add the respective definitions
SET( BUILD_EGFRvIII ON CACHE BOOL "Do you want to build EGFRvIII Estimation?" )
SET( BUILD_RECURRENCE ON CACHE BOOL "Do you want to build Recurrence Estimation?" )
SET( BUILD_SURVIVAL ON CACHE BOOL "Do you want to build Survival Estimation?" )
SET( BUILD_ISUBTYPE OFF CACHE BOOL "Do you want to build Imaging Subtype Estimation?" )
SET( BUILD_MSUBTYPE ON CACHE BOOL "Do you want to build Molecular Subtype Estimation?" )
SET( BUILD_ATLAS ON CACHE BOOL "Do you want to build Population Atlas?" )
SET( BUILD_GEODESIC ON CACHE BOOL "Do you want to build Geodesic Segmentation?" )
SET( BUILD_GEODESICTRAIN OFF CACHE BOOL "Do you want to build Geodesic Segmentation Train?" )
SET( BUILD_SBRT ON CACHE BOOL "Do you want to build SBRT Lung?" )
SET( BUILD_WHITESTRIPE ON CACHE BOOL "Do you want to build WhiteStripe Normalization?" )
SET( BUILD_ITKSNAP ON CACHE BOOL "Do you want to connect ITK-SNAP to CaPTk Console?" )
SET( BUILD_LIBRASINGLE ON CACHE BOOL "Do you want to connect LIBRA to CaPTk Console (needs MATLAB runtime)?" )
SET( BUILD_LIBRABATCH ON CACHE BOOL "Do you want to connect LIBRA to CaPTk Console (needs MATLAB runtime)?" )
SET( BUILD_CONFETTI ON CACHE BOOL "Do you want to connect Confetti to CaPTk Console (needs Python 2.7 and other dependencies)?" )
SET( BUILD_FETALBRAIN OFF CACHE BOOL "Do you want to build FetalBrain?" )
SET( BUILD_DEEPMEDIC OFF CACHE BOOL "Do you want to package DeepMedic?" )

IF(BUILD_EGFRvIII)
  ADD_DEFINITIONS(-DBUILD_EGFRvIII)
ENDIF()

IF(BUILD_RECURRENCE)
  ADD_DEFINITIONS(-DBUILD_RECURRENCE)
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

IF(BUILD_GEODESICTRAIN)
  ADD_DEFINITIONS(-DBUILD_GEODESICTRAIN)
ENDIF()

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src)
SET(INDIVIDUAL_APPLICATIONS_DIR ${PROJECT_SOURCE_DIR}/src/applications/individualApps)

# set stand-alone applications directory
SET( StandAloneCLIAppsFolder "CaPTk_CLIApps")

IF( BUILD_SBRT)
  SET( SBRT_PROJECT_NAME SBRT_Lung)
  SET ( Executable_1 ${SBRT_PROJECT_NAME}_Segment )
  SET ( Executable_2 ${SBRT_PROJECT_NAME}_Analyze )
  
  INCLUDE_DIRECTORIES( ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/segment/ )
  INCLUDE_DIRECTORIES( ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/analyze/ )
  INCLUDE_DIRECTORIES( ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/common_includes/)

  SET( segment_src 
    ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/segment/Joint_Segm.h
    ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/segment/Joint_Segm.cpp
  )

  SET( analyze_src 
    ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/analyze/LBPFeatures.h
    ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/analyze/LBPFeatures.cpp  
    ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/analyze/generateTextureFeatures.h
    ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/analyze/generateTextureFeatures.cpp
  )


  SET( Common_dependents
    #${CMAKE_CURRENT_SOURCE_DIR}/common_includes/cbicaPreProcessImage.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaLogging.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaLogging.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaUtilities.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaUtilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaCmdParser.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaCmdParser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaITKImageInfo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaITKImageInfo.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cbica_toolkit/cbicaITKSafeImageIO.h
  )
  
  ADD_EXECUTABLE(
    ${Executable_1} 
    ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/segment/main.cxx
    ${segment_src}
    ${Common_dependents}
  )
  #ADD_EXECUTABLE(
  #  ${Executable_2} 
  #  ${INDIVIDUAL_APPLICATIONS_DIR}/sbrt_source/src/analyze/main.cxx
  #  ${analyze_src}
  #  ${Common_dependents}
  #)

  SET_TARGET_PROPERTIES( ${Executable_1} PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )
  
  #SET_TARGET_PROPERTIES( ${Executable_2} PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )
  
  # Link the libraries to be used
  IF( ITKVtkGlue_LOADED )
    TARGET_LINK_LIBRARIES(
      ${Executable_1}
      ${Glue}  
      ${VTK_LIBRARIES} 
      ${ITK_LIBRARIES}
      ${OpenCV_LIBS}
    )
    #TARGET_LINK_LIBRARIES(
    #  ${Executable_2}
    #  ${Glue}  
    #  ${VTK_LIBRARIES} 
    #  ${ITK_LIBRARIES}
    #  ${OpenCV_LIBS}
    #)
  ELSE()
    TARGET_LINK_LIBRARIES(
      ${Executable_1}
      ${ITK_LIBRARIES}
      ${VTK_LIBRARIES} 
      ${OpenCV_LIBS}
    )
    #TARGET_LINK_LIBRARIES(
    #  ${Executable_2}
    #  ${ITK_LIBRARIES}
    #  ${VTK_LIBRARIES} 
    #  ${OpenCV_LIBS}
    #)
  ENDIF()

  INSTALL( TARGETS ${Executable_1} 
    #DESTINATION bin  
    BUNDLE DESTINATION .
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
  #INSTALL( TARGETS ${Executable_2} 
  #  #DESTINATION bin  
  #  BUNDLE DESTINATION .
  #  RUNTIME DESTINATION bin
  #  LIBRARY DESTINATION lib
  #  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  #  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  #)
ENDIF()

#IF( BUILD_WHITESTRIPE )
#  SET(WhiteStripeSrcs 
#    ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/main.cpp
#    ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/whiteStripe.cpp
#    ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/KernelFit.cpp)
#  SET(WhiteStripeHeaders
#    ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/testWhiteStripe.h
#    ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/whiteStripe.h
#    ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/KernelFit.h)
#
#  ADD_EXECUTABLE(WhiteStripe ${WhiteStripeSrcs} ${WhiteStripeHeaders})
#  
#  SET_TARGET_PROPERTIES( WhiteStripe PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )
#  
#  INSTALL( TARGETS WhiteStripe 
#      #DESTINATION bin  
#      BUNDLE DESTINATION .
#      RUNTIME DESTINATION bin
#      LIBRARY DESTINATION lib
#    CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
#    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
#  )
#  TARGET_LINK_LIBRARIES(WhiteStripe ${OpenCV_LIBS} ${ITK_LIBRARIES})
#
#  SET( BUILD_WHITESTRIPE_GUI ON CACHE BOOL "Do you want to build GUI?" )
#  IF( BUILD_WHITESTRIPE_GUI )
#    FIND_PACKAGE(Qt4 REQUIRED)
#    INCLUDE( ${QT_USE_FILE} )
#    SET(WhiteStripeGUISrc 
#      ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/whiteStripeGUI.cpp
#      ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/whiteStripe.cpp
#      ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/KernelFit.cpp
#      )
#    SET(WhiteStripeGUIIncludes 
#      ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/whiteStripe.h
#      ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/KernelFit.h
#      ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/utils.h
#      ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/whiteStripeGUI.h)
#    QT4_WRAP_CPP(MOCSrcs ${INDIVIDUAL_APPLICATIONS_DIR}/whitestripe_source/src/whiteStripeGUI.h)
#    
#    IF( NOT CAPTK_CLI_MODE )
#      ADD_EXECUTABLE(WhiteStripeGUI ${WhiteStripeGUISrc} ${WhiteStripeGUIIncludes} ${MOCSrcs})
#      
#      SET_TARGET_PROPERTIES( WhiteStripeGUI PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )
#  
#      TARGET_LINK_LIBRARIES(WhiteStripeGUI ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${OpenCV_LIBS})
#      INSTALL( TARGETS WhiteStripeGUI 
#        #DESTINATION bin  
#        BUNDLE DESTINATION .
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib
#        CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
#        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
#      )
#    ENDIF()
#    
#  ENDIF ()
#ENDIF()

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

    ADD_EXECUTABLE(Confetti ${ConfettiSrcs} ${ConfettiHeaders})
    TARGET_LINK_LIBRARIES(Confetti ${OpenCV_LIBS} ${ITK_LIBRARIES})
    SET_TARGET_PROPERTIES( Confetti PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )

    ADD_LIBRARY(ConfettiCore SHARED
      ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confetti.cpp
      ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confetti.h
      ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/multinomialModel.h
      ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/dtiUtils.h
      ${INDIVIDUAL_APPLICATIONS_DIR}/confetti_source/src/confettiApi.cpp
      )
    TARGET_LINK_LIBRARIES(ConfettiCore ${OpenCV_LIBS} ${ITK_LIBRARIES})

    # install all targets and scripts
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

ADD_DEFINITIONS(-DDATA_DIR="${DATA_DIR}")
ADD_DEFINITIONS(-DICON_DIR="${DATA_DIR}/icons")
ADD_DEFINITIONS(-DCAPTK_APP_LIST="${CAPTK_APP_LIST}")
ADD_DEFINITIONS(-DPREPROCESS_ALGOS="${PREPROCESSING_ALGOS}")
ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_GUI="${CAPTK_APP_LIST_PY_GUI}")
ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_CLI="${CAPTK_APP_LIST_PY_CLI}")

# generate moc files for each GUI element 
SET(GUI_SOURCES_DIR "${PROJECT_SOURCE_DIR}/src/gui")
SET(GUI_SOURCES_UPDATED ${GUI_SOURCES} )
FOREACH(gui_tool ${UI_ELEMENTS})
  
  IF( EXISTS "${GUI_SOURCES_DIR}/ui_${gui_tool}.h" )
    SET(
      GUI_SOURCES_UPDATED 
      ${GUI_SOURCES_UPDATED} 
      "${GUI_SOURCES_DIR}/ui_${gui_tool}.h"
      CACHE STRING "All UI elements" FORCE
    )
  ENDIF()
  IF( EXISTS "${GUI_SOURCES_DIR}/${gui_tool}.cpp" )
    SET(
      GUI_SOURCES_UPDATED 
      ${GUI_SOURCES_UPDATED} 
      "${GUI_SOURCES_DIR}/${gui_tool}.cpp"
      CACHE STRING "All UI elements" FORCE
    )
  ENDIF()
  IF( EXISTS "${GUI_SOURCES_DIR}/${gui_tool}.h" )
    SET(
      GUI_SOURCES_UPDATED 
      ${GUI_SOURCES_UPDATED} 
      "${GUI_SOURCES_DIR}/${gui_tool}.h"
      CACHE STRING "All UI elements" FORCE
    )
    
    QT4_WRAP_CPP(GUI_SOURCES_UPDATED 
    "${GUI_SOURCES_DIR}/${gui_tool}.h"
    )
  ENDIF()
  
ENDFOREACH(gui_tool)

# keep source files structured as folders for clarity
SET_PROPERTY( GLOBAL PROPERTY USE_FOLDERS ON )
SOURCE_GROUP( GUI FILES ${GUI_SOURCES_UPDATED} )
SOURCE_GROUP( CBICA-TK FILES ${CBICA-TK_SOURCES} )
SOURCE_GROUP( APPLICATIONS FILES ${APPLICATION_SOURCES} )
SOURCE_GROUP( APPLICATIONS_Dependencies FILES ${APPLICATION_DEPENDS} )

# set up the libraries for use in CaPTK
SET_PROPERTY( GLOBAL PROPERTY USE_FOLDERS ON )

SET( LIBNAME_GUI "${PROJECT_NAME}_GUI" )
SET( LIBNAME_Applications "${PROJECT_NAME}_Applications" )
SET( LIBNAME_FeatureExtractor "${PROJECT_NAME}_FeatureExtractor" )

# add all dependencies of ${application} to compilation
FILE( GLOB 
  FEATURE_EXTRACTION_DEPENDS
  "${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction_includes/*.h"
  "${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction_includes/*.cpp"
)

ADD_LIBRARY(
  ${LIBNAME_FeatureExtractor}
  ${FEATURE_EXTRACTION_DEPENDS}
  ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction.h
  ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction.cpp
  ${CBICA-TK_SOURCES}
)

TARGET_LINK_LIBRARIES(
  ${LIBNAME_FeatureExtractor}
  ${ITK_LIBRARIES}   
  ${OpenCV_LIBRARIES}
)      

ADD_LIBRARY( 
  ${LIBNAME_GUI} 
  ${GUI_SOURCES_UPDATED} 
)

ADD_LIBRARY( 
  ${LIBNAME_Applications} 
  ${APPLICATION_SOURCES}
  ${APPLICATION_DEPENDS}
  ${CBICA-TK_SOURCES}
)

TARGET_LINK_LIBRARIES(
  ${LIBNAME_Applications}
  ${ITK_LIBRARIES}   
  ${OpenCV_LIBRARIES}
)      

TARGET_LINK_LIBRARIES(
  ${LIBNAME_GUI}
  Qt4::QtCore
  Qt4::QtGui
  #QVTK         
  ${VTK_LIBRARIES}
  ${LIBNAME_FeatureExtractor}
)

SET(CaPTk_Libs_Folder "CaPTk_Libs")

ADD_DEPENDENCIES( ${LIBNAME_Applications} ${LIBNAME_GUI} ${LIBNAME_FeatureExtractor} )

SET_TARGET_PROPERTIES( ${LIBNAME_GUI} PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
SET_TARGET_PROPERTIES( ${LIBNAME_Applications} PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
SET_TARGET_PROPERTIES( ${LIBNAME_FeatureExtractor} PROPERTIES FOLDER "${CaPTk_Libs_Folder}" )
  
IF( NOT CAPTK_CLI_MODE )

  # define the executable
  IF(WIN32)
    # ensures no console pops up when launching the program
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP")
    ADD_EXECUTABLE(${EXE_NAME} WIN32 
      ${EXE_SOURCES}
      ${DATA_DIR}/icons/captk.rc
    )
  #ELSEIF(APPLE)
  #  # check how to get icon for mac
  #  ADD_EXECUTABLE(${EXE_NAME}
  #    ${EXE_SOURCES}
  #  )
  ELSE()
    ADD_EXECUTABLE(${EXE_NAME}
      ${EXE_SOURCES}
    )
  ENDIF()
  
  TARGET_LINK_LIBRARIES( ${EXE_NAME}
    ${DEPENDENT_LIBS} 
    ${LIBNAME_GUI}
    ${LIBNAME_Applications}
    ${LIBNAME_FeatureExtractor}
  )

  ADD_DEPENDENCIES( ${EXE_NAME} 
    ${LIBNAME_GUI} 
    ${LIBNAME_Applications}
    ${LIBNAME_FeatureExtractor}
  )
  
  # set properties for the executable
  SET_TARGET_PROPERTIES( ${EXE_NAME} 
    PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE
    OUTPUT_NAME ${EXE_NAME}
  )
  
ENDIF()

# find stand-alone CMD applications
FILE( GLOB STANDALONE_APPS
  "${PROJECT_SOURCE_DIR}/src/applications/*.cxx"
)
 
#MESSAGE( STATUS "[debug]STANDALONE_APPS == ${STANDALONE_APPS}")

FOREACH(application ${STANDALONE_APPS})

  GET_FILENAME_COMPONENT( STANDALONE_APP_NAME ${application} NAME_WE )

  ADD_EXECUTABLE(${STANDALONE_APP_NAME} 
    ${PROJECT_SOURCE_DIR}/src/applications/${STANDALONE_APP_NAME}.cxx 
    ${CBICA-TK_SOURCES}
    ${APPLICATION_DEPENDS}
  )
  
  TARGET_LINK_LIBRARIES(${STANDALONE_APP_NAME}
    ${DEPENDENT_LIBS} 
    ${LIBNAME_Applications}
    ${LIBNAME_FeatureExtractor}
  )

  SET_TARGET_PROPERTIES( ${STANDALONE_APP_NAME} PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )
  
  ADD_DEPENDENCIES( ${STANDALONE_APP_NAME} ${LIBNAME_GUI} ${LIBNAME_Applications} ) 
  
  INSTALL( TARGETS ${STANDALONE_APP_NAME} 
    BUNDLE DESTINATION .
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
ENDFOREACH()

IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/install" CACHE PATH "Default Install Path" FORCE )
ENDIF()

## Build test(s) only if CAPTK_PACKAGE_PROJECT is OFF
IF( NOT CAPTK_PACKAGE_PROJECT )
  SET( BUILD_TESTING FALSE CACHE BOOL "Build testing" )
  IF( BUILD_TESTING )
    SET (SAMPLE_DATA_FILE "CaPTk_SampleData.zip")
    IF (NOT EXISTS ${DATA_DIR}/README.txt) # check if data zip was extracted here before
      FILE( 
        DOWNLOAD "ftp://www.nitrc.org/var/home/groups/captk/downloads/CaPTk_SampleData.zip" 
        ${DATA_DIR}/${SAMPLE_DATA_FILE}
        SHOW_PROGRESS
        TIMEOUT 10000
      )
      EXECUTE_PROCESS(
        WORKING_DIRECTORY ${DATA_DIR}
        COMMAND ${CMAKE_COMMAND} -E tar xzf ${SAMPLE_DATA_FILE}
      )
      FILE( REMOVE ${DATA_DIR}/${SAMPLE_DATA_FILE} )
    ENDIF()
    ENABLE_TESTING()
    ADD_SUBDIRECTORY( testing )
    INCLUDE(CTest)
  ENDIF( BUILD_TESTING )
ENDIF()

# installation
IF (UNIX)
  IF(NOT LIBEXEC_INSTALL_DIR)
    SET(LIBEXEC_INSTALL_DIR "./libexec")
  ENDIF()
ENDIF()

INSTALL( TARGETS ${EXE_NAME} 
  BUNDLE DESTINATION .
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

## symlink creation for UNIX machines - currently only works for DEB or RPM 
#IF (UNIX OR APPLE)
  #INSTALL(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/bin/${EXE_NAME} ${CMAKE_INSTALL_PREFIX}/bin/captk)")
#ENDIF()

# ensure all licenses are in the package
FILE( GLOB_RECURSE ALL_LICENSES "${PROJECT_SOURCE_DIR}/licenses/*.txt" ) 
FOREACH( license ${ALL_LICENSES} )
  INSTALL( FILES "${license}" DESTINATION licenses )
ENDFOREACH()
#INSTALL(FILES "${README_FILE}" DESTINATION ".") # not needed any more

# ensure all icons are in the package
FILE( GLOB_RECURSE ALL_ICONS "${DATA_DIR}/icons/*.png" ) 
FOREACH( icon ${ALL_ICONS} )
  INSTALL(FILES "${icon}" DESTINATION data/icons )
ENDFOREACH()

SET( CAPTK_STYLESHEET_FILE "captk.qss" )
ADD_DEFINITIONS( -DCAPTK_STYLESHEET_FILE="${CAPTK_STYLESHEET_FILE}" )
#INSTALL(FILES "${DATA_DIR}/FetalModel.xml" DESTINATION data)
INSTALL(FILES "${DATA_DIR}/${CAPTK_STYLESHEET_FILE}" DESTINATION etc)
INSTALL(FILES "${DATA_DIR}/AAAC0_flair_pp_shrunk.nii.gz" DESTINATION data)
INSTALL(FILES "${DATA_DIR}/AAAC0_flair_pp_shrunk_testTumor.nii.gz" DESTINATION data)
INSTALL(DIRECTORY "${DATA_DIR}/sri24" DESTINATION data PATTERN ".svn" EXCLUDE)

IF(WIN32)
  FILE(GLOB_RECURSE ALL_DEPENDENT_LIBS "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/*.dll")
ELSE(WIN32)
  FILE(GLOB_RECURSE ALL_DEPENDENT_LIBS "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/*.so")
ENDIF()
FOREACH(dependent_lib ${ALL_DEPENDENT_LIBS})
  INSTALL(FILES "${dependent_lib}" DESTINATION bin)
ENDFOREACH()


# For new feature, add the feature definition name here; it needs to correspond to the file in "${DATA_DIR}/features"
SET( FeatureDefinitions "" )
FILE( GLOB allFeatureDefinitionFiles "${DATA_DIR}/features/*.csv" )
FOREACH(featureDefinitionFile ${allFeatureDefinitionFiles})

  GET_FILENAME_COMPONENT( CURRENT_DEF_FILE ${featureDefinitionFile} NAME_WE )

  INSTALL( 
    FILES "${featureDefinitionFile}" 
    DESTINATION data/features 
  )
  
  SET( FeatureDefinitions "${FeatureDefinitions} ${CURRENT_DEF_FILE}" )
ENDFOREACH()

INSTALL( FILES "${DATA_DIR}/features/batchmode/batch_featureextraction.csv" DESTINATION share/featureExtractionBatch )

SET( FeatureDefinitions ${FeatureDefinitions} CACHE STRING "Available Feature Definitions" FORCE)
ADD_DEFINITIONS(-DFeatureDefinitions="${FeatureDefinitions}")

# Install pre-trained models
#INSTALL(DIRECTORY "${DATA_DIR}/recurrence" DESTINATION data PATTERN ".svn" EXCLUDE)
#INSTALL(DIRECTORY "${DATA_DIR}/survival" DESTINATION data PATTERN ".svn" EXCLUDE)
#INSTALL(DIRECTORY "${PROJECT_SOURCE_DIR}/data/molecular" DESTINATION data PATTERN ".svn" EXCLUDE)
#INSTALL(DIRECTORY "${PROJECT_SOURCE_DIR}/data/subtype" DESTINATION data PATTERN ".svn" EXCLUDE)
#INSTALL(FILES "${DATA_DIR}/SBRT_SVM_Model.csv" DESTINATION data)

IF( WIN32 )
  INSTALL(DIRECTORY "${DATA_DIR}/GnuWin32" DESTINATION data PATTERN ".svn" EXCLUDE)
ENDIF()

IF(CAPTK_PACKAGE_PROJECT)

  SET(DEPENDENCY_DIRS ${VTK_LIBRARY_DIRS}/${CMAKE_CONFIGURATION_TYPES}/ ${QT_BINARY_DIR})
  IF(APPLE)
    SET(EXECUTABLE "${PROJECT_BINARY_DIR}/${CMAKE_CONFIGURATION_TYPES}/${EXE_NAME}.app")
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
  
  # [TBD] change this to a directory install, if possible 
  FILE(GLOB_RECURSE ALL_DOC_FILES "${PROJECT_SOURCE_DIR}/docs/tutorial/*.*")
  FOREACH(docFiles ${ALL_DOC_FILES})
    INSTALL(FILES "${docFiles}" DESTINATION share/doc)
  ENDFOREACH()  
  FILE(GLOB ALL_DOC_FILES "${PROJECT_BINARY_DIR}/docs/html/*.*")
  FOREACH(docFiles ${ALL_DOC_FILES})
    INSTALL(FILES "${docFiles}" DESTINATION share/doc)
  ENDFOREACH()
  FILE(GLOB ALL_DOC_FILES "${PROJECT_BINARY_DIR}/docs/html/search/*.*")
  FOREACH(docFiles ${ALL_DOC_FILES})
    INSTALL(FILES "${docFiles}" DESTINATION share/doc/search)
  ENDFOREACH()
  # ensure all icons are in the package
  FILE( GLOB ALL_IMGS "${PROJECT_SOURCE_DIR}/docs/images/gettingStarted/*.png" ) 
  FOREACH( img ${ALL_IMGS} )
    INSTALL(FILES "${img}" DESTINATION share/doc )
  ENDFOREACH()

  
  IF( WIN32 )
    INSTALL(FILES "${DATA_DIR}/icons/application/windows/captk.ico" DESTINATION share/icons)
  ELSEIF( UNIX )
    INSTALL( DIRECTORY "${DATA_DIR}/icons/application/linux/hicolor" DESTINATION share/icons)
  ENDIF()
  
  INSTALL(FILES "${DATA_DIR}/icons/captk.rc" DESTINATION data/icons)
  
  # [TBD] testing new script that checks system info
  INCLUDE(SystemSpecificInformation)
  MESSAGE( STATUS "SPECIFIC_COMPILER_NAME = ${SPECIFIC_COMPILER_NAME}" )
  MESSAGE( STATUS "SPECIFIC_SYSTEM_VERSION_NAME = ${SPECIFIC_SYSTEM_VERSION_NAME}" )
  MESSAGE( STATUS "SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR = ${SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR}" )
    
  # CPack settings
  IF(WIN32)
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_SOURCE_GENERATOR "ZIP")
    SET(CPACK_NSIS_DISPLAY_NAME "${PROJECT_NAME}") 
    SET(CPACK_NSIS_PACKAGE_NAME "${PROJECT_NAME}") 
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}_${${PROJECT_NAME}_Variant}\\\\${PROJECT_VERSION}")
    SET(CPACK_NSIS_INSTALL_ROOT "C:")
    SET(CPACK_NSIS_DEFINES "RequestExecutionLevel user")
    SET(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL "ON")
    SET(CPACK_NSIS_MUI_ICON "${DATA_DIR}/icons/application/windows/captk.ico")
    SET(CPACK_NSIS_MUI_FINISHPAGE_RUN "${EXE_NAME}.exe")
    #SET(CPACK_NSIS_INSTALLED_ICON_NAME "${PROJECT_NAME}${CMAKE_EXECUTABLE_SUFFIX}")
    #SET(CPACK_NSIS_MUI_FINISHPAGE_RUN "${CPACK_NSIS_INSTALLED_ICON_NAME}")
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
    SET(CPACK_NSIS_DELETE_ICONS "Delete '$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\${EXE_NAME}.lnk'")
    #SET(CPACK_NSIS_DELETE_ICONS "Delete '$DESKTOP\\\\${EXE_NAME}.lnk'")
    #SET(CPACK_NSIS_DELETE_ICONS_EXTRA "Delete '$INSTDIR\\\\${EXE_NAME}.lnk'")
    SET(CPACK_NSIS_DELETE_ICONS "Delete '$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\README.lnk'")
    SET(CPACK_NSIS_INSTALLED_ICON_NAME "${EXE_NAME}") 
    SET(CPACK_SYSTEM_NAME "Win64")
  ELSEIF(APPLE)
    SET(CMAKE_INSTALL_PREFIX "/Applications")
    SET(CPACK_GENERATOR "DragNDrop")
    SET(CPACK_DMG_FORMAT "UDBZ")
    SET(CPACK_DMG_VOLUME_NAME "${EXE_NAME}")
    SET(CPACK_SYSTEM_NAME "OSX")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}_${${PROJECT_NAME}_Variant}/${PROJECT_VERSION}")
    #SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
    #SET(CPACK_DMG_DS_STORE "${ICONS_DIR}/DMGDSStore")
    #SET(CPACK_DMG_BACKGROUND_IMAGE "${ICONS_DIR}/DMGBackground.png")
  ELSEIF(UNIX)
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}_${${PROJECT_NAME}_Variant}/${PROJECT_VERSION}")
    #SET(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
    SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
    SET(CPACK_GENERATOR "STGZ;DEB")
    
    ## if LIBRA ever moves to C++, use this for packaging for Debian
    #IF(${SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR} STREQUAL "DEB")
    #  SET(CPACK_GENERATOR "STGZ;DEB")
    #  SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "${PROJECT_VENDOR}")
    #ENDIF()
    
    SET(CPACK_SOURCE_GENERATOR "TGZ")
    SET(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
    SET(CPACK_SYSTEM_NAME "Linux-i686")
    SET(CPACK_TOPLEVEL_TAG "Linux-i686")
    SET(CPACK_STRIP_FILES "${EXE_NAME}")
  ENDIF()
  # common
  SET(CPACK_SOURCE_STRIP_FILES TRUE )
  SET(CPACK_PACKAGE_EXECUTABLES "${EXE_NAME}" "${EXE_NAME}")
  SET(CPACK_CREATE_DESKTOP_LINKS "${EXE_NAME}" ) 
  SET(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}_${PROJECT_VERSION}")
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
  set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/bin/;{CMAKE_SOURCE_DIR}/build/;${CMAKE_SOURCE_DIR}/.svn/")

  INCLUDE(InstallRequiredSystemLibraries)
  INCLUDE(CPack)

  IF(EXISTS ${EXE_NAME}) # don't do anything if the EXE is not found - doing this instead of "install_qt4_exe" since it works
    # this will no longer be needed with qt5 since it comes with its own packaging tool
    #IF(APPLE OR WIN32)
      #INCLUDE(DeployQt4)
      #INCLUDE(BundleUtilities)
      #MESSAGE(STATUS "Qt libs: ${QT_QTCORE_LIBRARY};${QT_QTGUI_LIBRARY}")
      #MESSAGE(STATUS "Qt lib dir: ${QT_LIBRARY_DIR}")
      #MESSAGE(STATUS "Qt bin dir: ${QT_BINARY_DIR}")
      #MESSAGE(STATUS "VTK bin dir: ${VTK_LIBRARY_DIRS}/${CMAKE_CONFIGURATION_TYPES}/")
      #FIXUP_BUNDLE("${EXECUTABLE}" "" "${DEPENDENCY_DIRS}")
      #FIXUP_QT4_EXECUTABLE("${EXECUTABLE}" "" "" "${QT_LIBRARY_DIR}" "" "")
      #SET(qtconf_dest_dir /Contents/Resources)
      #INSTALL(CODE 
      #  "
      #  FILE(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"\")
      #  " 
      #  COMPONENT Runtime
      #)
      IF(WIN32)
        INSTALL( CODE 
          "
          INCLUDE(BundleUtilities)
          FIXUP_BUNDLE( \"${EXE_NAME}\" \"\" \"${QT_BINARY_DIR}\" )
          " 
          Component Runtime 
        )
      #INSTALL_QT4_EXECUTABLE("${EXECUTABLE}" "" "${QT_LIBS}" "${QT_LIBRARY_DIR}" "" "" "")
      ENDIF()
    #ENDIF()
  ENDIF()

ENDIF(CAPTK_PACKAGE_PROJECT)

ENDIF()