PROJECT(CaPTk_GUI)

FIND_PACKAGE(Qt5 COMPONENTS Core Widgets Concurrent Test REQUIRED)

SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY 
	${CaPTk_BINARY_DIR})

ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_GUI="${CAPTK_APP_LIST_PY_GUI}")
ADD_DEFINITIONS(-DCAPTK_APP_LIST_PY_CLI="${CAPTK_APP_LIST_PY_CLI}")
ADD_DEFINITIONS(-DFeatureDefinitions="${FeatureDefinitions}")

SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES( 	${CMAKE_CURRENT_SOURCE_DIR}/inc
                      ${CaPTk_SOURCE_DIR}/src/applications/FeatureExtraction/src
                     )
  
SET( ${PROJECT_NAME}_HDRS
  inc/ApplicationBase.h
  inc/ASyncAppManager.h
  inc/fBottomImageInfoTip.h
  inc/fDCM2NIfTI.h
  inc/fDeepMedicDialog.h
  inc/fDeepMedicNormDialog.h
  inc/fDiffusionMeasuresDialog.h
  inc/fDirectionalityDialog.h
  inc/fDrawingPanel.h
  inc/fEGFRvIIIDialog.h
  inc/fFeatureDialog.h
  inc/fFeaturePanel.h
  inc/fFetalBrain.h
  inc/fHelpAppDialog.h
  inc/fHelpDialog.h
  inc/fHelpTutorial.h
  inc/fHistoMatchDialog.h
  inc/fImagesPanel.h
  inc/fImagingSubtypeDialog.h
  inc/fLungFieldDialog.h
  inc/fMainWindow.h
  inc/fMolecularSubtypeDialog.h
  inc/fPCADialog.h
  inc/fPerfusionMeasuresDialog.h
  inc/fPopulationAtlasDialog.h
  inc/fPreprocessingDialog.h
  inc/fProgressDialog.h
  inc/fPseudoProgressionDialog.h
  inc/fRecurrenceDialog.h
  inc/fRegistrationDialog.h
  inc/fSBRTAnalysisDialog.h
  inc/fSBRTNoduleDialog.h
  inc/fSegmentationPanel.h
  inc/fSkullStripDialog.h
  inc/fSurvivalDialog.h
  inc/fTrainingDialog.h
  inc/fTumorPanel.h
  inc/fWhiteStripeDialog.h
  inc/InteractorStyleNavigator.h
  inc/Landmarks.h
  inc/OutputInteractorStyleNavigator.h
  inc/QBorderWidget.h
  inc/QTablePushButton.h
  inc/SimpleImageManager.h
  inc/Slicer.h
  inc/SlicerManager.h
  inc/SlicerManagerCommand.h
  inc/ui_ApplicationBase.h
  inc/ui_fBottomImageInfoTip.h
  inc/ui_fDCM2NIfTI.h
  inc/ui_fDeepMedicDialog.h
  inc/ui_fDeepMedicNormDialog.h
  inc/ui_fDiffusionMeasuresDialog.h
  inc/ui_fDirectionalityDialog.h
  inc/ui_fDrawingPanel.h
  inc/ui_fEGFRvIIIDialog.h
  inc/ui_fFeatureDialog.h
  inc/ui_fFeaturePanel.h
  inc/ui_fFetalBrain.h
  inc/ui_fHelpAppDialog.h
  inc/ui_fHelpDialog.h
  inc/ui_fHistoMatchDialog.h
  inc/ui_fImagesPanel.h
  inc/ui_fImagingSubtypeDialog.h
  inc/ui_fLungFieldDialog.h
  inc/ui_fMainWindow.h
  inc/ui_fMolecularSubtypeDialog.h
  inc/ui_fPCADialog.h
  inc/ui_fPerfusionMeasuresDialog.h
  inc/ui_fPopulationAtlasDialog.h
  inc/ui_fPreprocessingDialog.h
  inc/ui_fProgressDialog.h
  inc/ui_fPseudoProgressionDialog.h
  inc/ui_fRecurrenceDialog.h
  inc/ui_fRegistrationDialog.h
  inc/ui_fSBRTAnalysisDialog.h
  inc/ui_fSBRTNoduleDialog.h
  inc/ui_fSegmentationPanel.h
  inc/ui_fSkullStripDialog.h
  inc/ui_fSurvivalDialog.h
  inc/ui_fTrainingDialog.h
  inc/ui_fTumorPanel.h
  inc/ui_fWhiteStripeDialog.h
  inc/ui_QTablePushButton.h
  inc/ui_SlicerManager.h
  "${CaPTk_SOURCE_DIR}/src/applications/GeodesicSegmentation.h"
)

SET( ${PROJECT_NAME}_SRCS
  src/ApplicationBase.cpp
  src/ASyncAppManager.cpp
  src/fBottomImageInfoTip.cpp
  src/fDCM2NIfTI.cpp
  src/fDeepMedicDialog.cpp
  src/fDeepMedicNormDialog.cpp
  src/fDiffusionMeasuresDialog.cpp
  src/fDirectionalityDialog.cpp
  src/fDrawingPanel.cpp
  src/fEGFRvIIIDialog.cpp
  src/fFeatureDialog.cpp
  src/fFeaturePanel.cpp
  src/fFetalBrain.cpp
  src/fHelpAppDialog.cpp
  src/fHelpDialog.cpp
  src/fHelpTutorial.cpp
  src/fHistoMatchDialog.cpp
  src/fImagesPanel.cpp
  src/fImagingSubtypeDialog.cpp
  src/fLungFieldDialog.cpp
  src/fMainWindow.cpp
  src/fMolecularSubtypeDialog.cpp
  src/fPCADialog.cpp
  src/fPerfusionMeasuresDialog.cpp
  src/fPopulationAtlasDialog.cpp
  src/fPreprocessingDialog.cpp
  src/fProgressDialog.cpp
  src/fPseudoProgressionDialog.cpp
  src/fRecurrenceDialog.cpp
  src/fRegistrationDialog.cpp
  src/fSBRTAnalysisDialog.cpp
  src/fSBRTNoduleDialog.cpp
  src/fSegmentationPanel.cpp
  src/fSkullStripDialog.cpp
  src/fSurvivalDialog.cpp
  src/fTrainingDialog.cpp
  src/fTumorPanel.cpp
  src/fWhiteStripeDialog.cpp
  src/InteractorStyleNavigator.cpp
  src/Landmarks.cpp
  src/QTablePushButton.cpp
  src/SimpleImageManager.cpp
  src/Slicer.cpp
  src/SlicerManager.cpp
  src/SlicerManagerCommand.cpp
  "${CaPTk_SOURCE_DIR}/src/applications/GeodesicSegmentation.cpp"
)

QT5_WRAP_CPP(${PROJECT_NAME}_MOC_HDRS
	${${PROJECT_NAME}_HDRS}
)

INCLUDE_DIRECTORIES( ${CaPTk_SOURCE_DIR}/src
                     ${CaPTk_SOURCE_DIR}/src/cbica_toolkit
                     ${CaPTk_SOURCE_DIR}/src/applications/common_includes
                     ${CaPTk_SOURCE_DIR}/src/applications
                     ${CaPTk_SOURCE_DIR}/src/applications/FeatureExtraction_includes
                     ${CaPTk_SOURCE_DIR}/src/applications/WhiteStripe_includes
                     ${CaPTk_SOURCE_DIR}/src/Common/inc
                     ${CaPTk_SOURCE_DIR}/src/applications/FeatureExtraction/src/depends
                     )
                     
ADD_LIBRARY(
  ${PROJECT_NAME}
  ${${PROJECT_NAME}_HDRS}
  ${${PROJECT_NAME}_SRCS}
  ${${PROJECT_NAME}_MOC_HDRS}
)

TARGET_LINK_LIBRARIES(
  ${PROJECT_NAME}
  ${DEPENDENT_LIBS}
  ${LIBNAME_CBICATK}
  ${LIBNAME_FeatureExtractor}
)