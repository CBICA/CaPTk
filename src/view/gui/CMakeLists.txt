# This CMakeLists file will add the source files in this directory for compilation

#find_package(Qt4 REQUIRED)

PROJECT(CaPTk_GUI)

# NEW CHANGES
FIND_PACKAGE(Qt5 COMPONENTS Core Widgets Concurrent Test REQUIRED)

SET(CMAKE_AUTOMOC ON)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES( ${CACHED_INCLUDE_DIRS} )

SET( GUI
  ApplicationBase
  SlicerManager
  QTablePushButton
  fMainWindow
  fProgressDialog
  fBottomImageInfoTip
  fTumorPanel
  fImagesPanel
  fDrawingPanel
  fFeaturePanel
  fFeatureDialog
  fRecurrenceDialog
  fSurvivalDialog
  fPseudoProgressionDialog
  fPopulationAtlasDialog
  fMolecularSubtypeDialog
  fEGFRvIIIDialog
  fImagingSubtypeDialog
  fFetalBrain
  fSkullStripDialog
  fPCADialog
  fTrainingDialog
  fPerfusionMeasuresDialog
  fPerfusionAlignmentDialog
  fDiffusionMeasuresDialog
  fDCM2NIfTI
  fHistoMatchDialog
  fWhiteStripeDialog
  fRegistrationDialog
  fPreprocessingDialog
  fHelpDialog
  fHelpTutorial
  fDirectionalityDialog
  fDeepMedicDialog
  fDeepMedicNormDialog
  fSBRTNoduleDialog
  fSBRTAnalysisDialog
  fTexturePipelineDialog
  AppearancePage
  CaPTkDockWidget
  fBiasCorrectionDialog
  # fAppDownloadDialog
  fBraTSSegmentation
  #fSegmentationPanel
  #fHelpAppDialog
)


SET( GUI_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/SimpleImageManager.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/SimpleImageManager.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Slicer.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Slicer.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/SlicerManagerCommand.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/SlicerManagerCommand.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/ComparisonViewerCommand.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/ComparisonViewerCommand.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Landmarks.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Landmarks.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/InteractorStyleNavigator.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/InteractorStyleNavigator.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/QBorderWidget.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/OutputInteractorStyleNavigator.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/ASyncAppManager.h"
  # "${CMAKE_CURRENT_SOURCE_DIR}/ThreadedExtraction.h"
  # "${CMAKE_CURRENT_SOURCE_DIR}/ThreadedExtraction.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/CheckOpenGLVersion.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/CheckOpenGLVersion.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/PreferencesDialog.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/PreferencesDialog.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/IPreferencePage.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/IPreferencePage.cpp"
  # "${CMAKE_CURRENT_SOURCE_DIR}/ApplicationPreferences.h"
  # "${CMAKE_CURRENT_SOURCE_DIR}/ApplicationPreferences.cpp"
  # "${CMAKE_CURRENT_SOURCE_DIR}/QZipReader.h"
  # "${CMAKE_CURRENT_SOURCE_DIR}/QZipWriter.h"
  # "${CMAKE_CURRENT_SOURCE_DIR}/QZip.cpp"
  # "${CMAKE_CURRENT_SOURCE_DIR}/ApplicationDownloadManager.h"
  # "${CMAKE_CURRENT_SOURCE_DIR}/ApplicationDownloadManager.cpp"
  "${CaPTk_SOURCE_DIR}/src/applications/GeodesicSegmentation.h"
  "${CaPTk_SOURCE_DIR}/src/applications/GeodesicSegmentation.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/SystemInformationDisplayWidget.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/SystemInformationDisplayWidget.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/SystemInformation.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/SystemInformation.cpp"
  #"${CMAKE_CURRENT_SOURCE_DIR}/ApplicationBase.h"
  #"${CMAKE_CURRENT_SOURCE_DIR}/NearFarDrawingManager.h"
  #"${CMAKE_CURRENT_SOURCE_DIR}/NearFarDrawingManager.cpp"
)

# save file list to cache
SET(GUI_SOURCES ${GUI_SOURCES}
  CACHE STRING "GUI sources")

SET( UI_ELEMENTS ${GUI}
  CACHE STRING "UI Elements"
)

MARK_AS_ADVANCED(${GUI_SOURCES})
MARK_AS_ADVANCED(${UI_ELEMENTS})

SET( CAPTK_STYLESHEET_FILE "captk.qss" )
ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_GUI="${CAPTK_APP_LIST_PY_GUI}")
ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_CLI="${CAPTK_APP_LIST_PY_CLI}")
ADD_DEFINITIONS(-DFeatureDefinitions="${FeatureDefinitions}")
ADD_DEFINITIONS( -DCAPTK_STYLESHEET_FILE="${CAPTK_STYLESHEET_FILE}" )

#SET(GUI_SOURCES_DIR "${PROJECT_SOURCE_DIR}/src/gui")
SET(GUI_SOURCES_DIR "${PROJECT_SOURCE_DIR}")
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

    QT5_WRAP_CPP(GUI_SOURCES_UPDATED
    "${GUI_SOURCES_DIR}/${gui_tool}.h"
    )
  ENDIF()

ENDFOREACH(gui_tool)

INCLUDE_DIRECTORIES( ${CaPTk_SOURCE_DIR}/src
  ${CaPTk_SOURCE_DIR}/src/cbica_toolkit/src
  ${CaPTk_SOURCE_DIR}/src/applications/common_includes
  ${CaPTk_SOURCE_DIR}/src/applications
  ${CaPTk_SOURCE_DIR}/src/applications/FeatureExtraction_includes
  ${CaPTk_SOURCE_DIR}/src/applications/WhiteStripe_includes
  ${CaPTk_SOURCE_DIR}/src/common/inc
  ${CaPTk_SOURCE_DIR}/src/applications/FeatureExtraction/src/depends
  ${CaPTk_SOURCE_DIR}/src/cbica_toolkit/src/IO/inc
  ${CaPTk_SOURCE_DIR}/src/applications/Preprocessing/src/depends
  #${CMAKE_CURRENT_SOURCE_DIR}/inc
  #${VTK_INCLUDE_DIRS} 
)
                     
ADD_LIBRARY(
  #${LIBNAME_GUI}
  ${PROJECT_NAME}
  ${GUI_SOURCES_UPDATED}
)

TARGET_LINK_LIBRARIES(
  ${PROJECT_NAME}
  #${LIBNAME_GUI}
  GeodesicTrainingCaPTkLib
  ${DEPENDENT_LIBS}
  ${LIBNAME_CBICATK}
  ${LIBNAME_FeatureExtractor}
  #CaPTk_IO
)
