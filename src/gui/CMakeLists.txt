# This CMakeLists file will add the source files in this directory for compilation

find_package(Qt4 REQUIRED)
SET(CMAKE_AUTOMOC ON)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

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
  fPopulationAtlasDialog
  fMolecularSubtypeDialog
  fImagingSubtypeDialog
  fFetalBrain
  fSkullStripDialog
  fPCADialog
  fPerfusionMeasuresDialog
  fDiffusionMeasuresDialog
  fPerfusionAlignmentDialog
  fDCM2NIfTI
  fHistoMatchDialog
  fWhiteStripeDialog
  fRegistrationDialog
  fPreprocessingDialog
  fHelpDialog
  fHelpTutorial
  fDirectionalityDialog
  #fHelpAppDialog  
)

SET( GUI_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/SimpleImageManager.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/SimpleImageManager.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Slicer.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Slicer.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/SlicerManagerCommand.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/SlicerManagerCommand.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/Landmarks.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/Landmarks.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/InteractorStyleNavigator.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/InteractorStyleNavigator.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/QBorderWidget.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/OutputInteractorStyleNavigator.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/ASyncAppManager.h"
  "${PROJECT_SOURCE_DIR}/src/applications/GeodesicSegmentation.h"
  "${PROJECT_SOURCE_DIR}/src/applications/GeodesicSegmentation.cpp"
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