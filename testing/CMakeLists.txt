# add source code include directories
INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/src/
  ${PROJECT_SOURCE_DIR}/src/applications/
  ${PROJECT_SOURCE_DIR}/src/applications/py
  ${PROJECT_SOURCE_DIR}/src/applications/common_includes/
  ${PROJECT_SOURCE_DIR}/src/cbica_toolkit/
  ${APPLICATION_INCLUDES}
)

SET( TEST_EXE_NAME CaPTk_Tests )

SET( TESTING_DATA_DIR ${CMAKE_CURRENT_BINARY_DIR}/TestData )
ADD_DEFINITIONS( -DTEST_DATA_DIR=${TESTING_DATA_DIR} )

IF( NOT EXISTS "${TESTING_DATA_DIR}" )  
  
  SET( LFS_FILE_TO_CHECK "${PROJECT_SOURCE_DIR}/data/TestData.zip" )
  SET( FILE_TO_EXTRACT "${CMAKE_CURRENT_BINARY_DIR}/TestData.zip")
  SET( DOWNLOAD_URL "ftp://www.nitrc.org/home/groups/captk/downloads/SampleData/TestData.zip" ) 
  
  IF( NOT EXISTS "${FILE_TO_EXTRACT}" )
    
    # copy from LFS folder
    IF( EXISTS ${LFS_FILE_TO_CHECK} )
      CONFIGURE_FILE( ${LFS_FILE_TO_CHECK} ${FILE_TO_EXTRACT} )
    ENDIF()

    # do not re-download if the LFS fetch worked
    IF(NOT EXISTS ${FILE_TO_EXTRACT})
      MESSAGE( STATUS "Downloading all the sample data using following link: ${DOWNLOAD_URL}" )
      FILE(DOWNLOAD "${DOWNLOAD_URL}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
      IF(NOT STATUS_CODE EQUAL 0)
        MESSAGE(FATAL_ERROR "Failed to download TestData. Status=${STATUS_CODE}")
      ENDIF()
    ENDIF()
    
  ENDIF()
  
  FILE(MAKE_DIRECTORY ${TESTING_DATA_DIR})
  MESSAGE( STATUS "Extracting all sample data" )
  
  IF( EXISTS "${FILE_TO_EXTRACT}" )

    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      RESULT_VARIABLE RESULT_CODE
    )
    
    IF(NOT RESULT_CODE EQUAL 0)
      MESSAGE( WARNING "Extracting the test data failed" )
    ENDIF()
    
  ENDIF()
  
ENDIF()

SET( TESTING_OUTPUT_DIR ${TESTING_DATA_DIR}/outputs )
FILE(MAKE_DIRECTORY ${TESTING_OUTPUT_DIR})

SET(CMAKE_CXX_STANDARD 11 CACHE STRING "C++ standard to be used" )

#set(CMAKE_AUTOMOC ON)    
 
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
    
    QT5_WRAP_CPP(GUI_SOURCES_UPDATED 
    "${GUI_SOURCES_DIR}/${gui_tool}.h"
    )
  ENDIF()
  
ENDFOREACH(gui_tool)

ADD_EXECUTABLE(${TEST_EXE_NAME} AllTests.cxx )

FIND_PACKAGE(Qt5 COMPONENTS Core Gui Svg Widgets WebView WebEngine WebEngineCore WebEngineWidgets Concurrent Test REQUIRED)
LINK_DIRECTORIES(${QT_LIBRARY_DIR})
#INCLUDE(${QT_USE_FILE})
SET(QT_USE_QTGUI TRUE)
SET(QT_USE_QTMAIN TRUE) 
SET(QT_USE_QTGUI TRUE)
ADD_DEFINITIONS(${QT_DEFINITIONS})
ADD_DEFINITIONS(-DQT_PLUGIN)
ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)
LINK_DIRECTORIES(${QT_LIBRARY_DIR})

TARGET_LINK_LIBRARIES(${TEST_EXE_NAME} ${DEPENDENT_LIBS} ${LIBNAME_CBICATK} ${LIBNAME_FeatureExtractor} ${LIBNAME_GUI} ${LIBNAME_Applications} )

#ADD_TEST(NAME TestNumber1 COMMAND CaPTk_Console --run_test "random")

###
# Testing Suite
###

# NOTE: CMake treats any output to stdout as a failure (why???) so testing outputs such as -h -v or logging options simply result in
# failure since printing a usage to stdout is considered a failed state.

# Also note that COMMAND <program> expects the program to be in the path so for any type of, say, CI
# make sure to export the build path before running tests or you'll sit there wondering why 100% of tests fail

# IN ADDITION. TEST_EXE_NAME refers to ALLTESTS.cxx so, pretty much anything that uses that will fail

# DeepMedic Tests

# std::system((getApplicationPath("DeepMedic") +
# " -t1" + file_t1Norm +
# " -t1c " + file_t1ceNorm +
# " -t2 " + file_t2Norm +
# " -fl " + file_flNorm +
# " -o " + outputDirectory).c_str());

# DeepMedic Tests
# ADD_TEST(NAME DeepMedicTestSeg COMMAND DeepMedic -md ${PROJECT_SOURCE_DIR}/data/deepMedic/saved_models/brainTumorSegmentation -o ${TESTING_OUTPUT_DIR}/deepmedic/ -t1 ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_t1_LPS_r.nii.gz -t1c ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_t1ce_LPS_r.nii.gz -t2 ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_t2_LPS_r.nii.gz -fl ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_flair_LPS_r.nii.gz )
# ADD_TEST(NAME DeepMedicTestSkull COMMAND DeepMedic -md ${PROJECT_SOURCE_DIR}/data/deepMedic/saved_models/skullStripping -o ${TESTING_OUTPUT_DIR}/deepmedic/ -t1 ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_t1_LPS_r.nii.gz -t1c ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_t1ce_LPS_r.nii.gz -t2 ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_t2_LPS_r.nii.gz -fl ${TESTING_DATA_DIR}/deepmedic/AAAB_2006.10.28_flair_LPS_r.nii.gz )

# Diffusion Derivatives (Incomplete)

ADD_TEST(NAME DiffusionDerivativesHelpTest COMMAND DiffusionDerivatives -h )
ADD_TEST(NAME DiffusionDerivativesVersionTest COMMAND DiffusionDerivatives -v )

# Directionality Estimate Tests

# ADD_TEST(NAME DirectionalityEstimateTest COMMAND DirectionalityEstimate -l ${TESTING_DATA_DIR}/direst/img_timepoint_1.nii.gz -o ${TESTING_OUTPUT_DIR}/dirEst.txt )

ADD_TEST(NAME DirectionalityEstimateHelpTest COMMAND DirectionalityEstimate -h )
ADD_TEST(NAME DirectionalityEstimateVersionTest COMMAND DirectionalityEstimate -v )

# EGFRvIIIIndexPredictor Tests

ADD_TEST(NAME EGFRvIIIIndexPredictorNewModelTest COMMAND EGFRvIIIIndexPredictor -i ${TESTING_DATA_DIR}/EGFRvIIISVMIndex_1.6.1 -t 0 -o ${TESTING_OUTPUT_DIR}/EGFRvIIIIndexPredictor )
ADD_TEST(NAME EGFRvIIIIndexPredictorModelTest COMMAND EGFRvIIIIndexPredictor -i ${TESTING_DATA_DIR}/EGFRvIIISVMIndex_1.6.1 -t 1 -o ${TESTING_OUTPUT_DIR}/EGFRvIIIIndexPredictor -m ${TESTING_OUTPUT_DIR}/EGFRvIIIIndexPredictor)

ADD_TEST(NAME EGFRvIIIIndexPredictorHelpTest COMMAND EGFRvIIIIndexPredictor -h )
ADD_TEST(NAME EGFRvIIIIndexPredictorVersionTest COMMAND EGFRvIIIIndexPredictor -v )

# EGFRvIIISurrogateIndex Tests

ADD_TEST(NAME EGFRvIIISurrogateDSCMRITest COMMAND EGFRvIIISurrogateIndex -i ${TESTING_DATA_DIR}/egfrsur/DSC-MRI_data.nii.gz -m ${TESTING_DATA_DIR}/egfrsur/Near_Far_masks.nii.gz )
ADD_TEST(NAME EGFRvIIISurrogateT1GDTest COMMAND EGFRvIIISurrogateIndex -i ${TESTING_DATA_DIR}/egfrsur/T1-Gd_data.nii.gz -m ${TESTING_DATA_DIR}/egfrsur/Near_Far_masks.nii.gz )
ADD_TEST(NAME EGFRvIIISurrogateT2FlairTest COMMAND EGFRvIIISurrogateIndex -i ${TESTING_DATA_DIR}/egfrsur/T2-FLAIR_data.nii.gz -m ${TESTING_DATA_DIR}/egfrsur/Near_Far_masks.nii.gz )

ADD_TEST(NAME EGFRvIIISurrogateIndexHelpTest COMMAND EGFRvIIISurrogateIndex -h )
ADD_TEST(NAME EGFRvIIISurrogateIndexVersionTest COMMAND EGFRvIIISurrogateIndex -v )

# Feature Extraction Tests
# ADD_TEST(NAME FeatureExtractionRegressionTest COMMAND FeatureExtraction -n Test -i ${TESTING_DATA_DIR}/FeatureExtraction/image.nii.gz -p ${TESTING_DATA_DIR}/FeatureExtraction/params.csv -m ${TESTING_DATA_DIR}/FeatureExtraction/mask.nii.gz -o ${TESTING_OUTPUT_DIR}/featExMapWrite_ref.csv -t Mod -r 1 -l ROI --verticalConc 1 -ut ${TESTING_DATA_DIR}/FeatureExtraction/baseline.csv )
ADD_TEST(NAME FeatureExtractionMultipleROITest COMMAND FeatureExtraction -n Test -i ${TESTING_DATA_DIR}/FeatureExtraction/image.nii.gz -p ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/data/1_params_default.csv -m ${TESTING_DATA_DIR}/FeatureExtraction/mask.nii.gz -o ${TESTING_OUTPUT_DIR}/featExMapWrite_lattice.csv -t Mod -r 1,2 -l ED,NC )
ADD_TEST(NAME FeatureExtractionMultipleModalityTest COMMAND FeatureExtraction -n Test -i ${TESTING_DATA_DIR}/FeatureExtraction/image.nii.gz,${TESTING_DATA_DIR}/FeatureExtraction/image.nii.gz -p ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/data/1_params_default.csv -m ${TESTING_DATA_DIR}/FeatureExtraction/mask.nii.gz -o ${TESTING_OUTPUT_DIR}/featExMapWrite_multiMod.csv -t Mod1,Mod2 -r 1 -l ED )
ADD_TEST(NAME FeatureExtractionMultipleModalityMultipleROITest COMMAND FeatureExtraction -n Test -i ${TESTING_DATA_DIR}/FeatureExtraction/image.nii.gz,${TESTING_DATA_DIR}/FeatureExtraction/image.nii.gz -p ${PROJECT_SOURCE_DIR}/src/applications/FeatureExtraction/data/1_params_default.csv -m ${TESTING_DATA_DIR}/FeatureExtraction/mask.nii.gz -o ${TESTING_OUTPUT_DIR}/featExMapWrite_multiModMultiROI.csv -t Mod1,Mod2 -r 1,2 -l ED,NC )
# ADD_TEST(NAME FeatureExtractionBatchTest COMMAND FeatureExtraction )

ADD_TEST(NAME FeatureExtractionHelpTest COMMAND FeatureExtraction -h )
ADD_TEST(NAME FeatureExtractionVersionTest COMMAND FeatureExtraction -v )

# Geodesic Segmentation Tests (Incomplete)

ADD_TEST(NAME GeodesicSegmentationHelpTest COMMAND GeodesicSegmentation -h )
ADD_TEST(NAME GeodesicSegmentationVersionTest COMMAND GeodesicSegmentation -v )

# Greedy Registration Tests (Incomplete)

# -reg -trf -i moving.nii.gz -f fixed.nii.gz -o output.nii.gz -t matrix.mat -a -m MI -n 100x50x5 -th 4

ADD_TEST(NAME GreedyRegistrationHelpTest COMMAND GreedyRegistration -h )
ADD_TEST(NAME GreedyRegistrationVersionTest COMMAND GreedyRegistration -v )

# Imaging Subtype Predictor Tests

# Will not be in the latest 2018 RSNA release, not being tested
# ADD_TEST(NAME ImagingSubtypePredictorTest COMMAND ImagingSubtypePredictor -i ${TESTING_DATA_DIR}/subtype -o ${TESTING_OUTPUT_DIR}/subtype/ )

# ADD_TEST(NAME ImagingSubtypeHelpTest COMMAND ImagingSubtype -h )
# ADD_TEST(NAME ImagingSubtypeVersionTest COMMAND ImagingSubtype -v )

# Molecular Subtype Predictor Tests

ADD_TEST(NAME MolecularSubtypePredictorT0Test COMMAND MolecularSubtypePredictor -t 0 -i ${TESTING_DATA_DIR}/molsub/test_samples/AAUA -o ${TESTING_OUTPUT_DIR}/molsub/ )
ADD_TEST(NAME MolecularSubtypePredictorT1Test COMMAND MolecularSubtypePredictor -t 1 -i ${TESTING_DATA_DIR}/molsub/test_samples/AAUA -o ${TESTING_OUTPUT_DIR}/molsub/ -m ${TESTING_OUTPUT_DIR}/molsub/)

ADD_TEST(NAME MolecularSubtypePredictorHelpTest COMMAND MolecularSubtypePredictor -h )
ADD_TEST(NAME MolecularSubtypePredictorVersionTest COMMAND MolecularSubtypePredictor -v )

# PerfusionDerivatives Tests

# NOTE: I am very bothered by the fact that PerfusionDerivatives -r and -p flags be not only required arguments
# but apparently they *also* have default values which completely defeats the purpose of them being required.

ADD_TEST(NAME PerfusionDerivativesP1R1PH1Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 1 -o ${TESTING_OUTPUT_DIR}/perf -p 1 -r 1 -pH 1)
ADD_TEST(NAME PerfusionDerivativesP1R1PH0Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 1 -o ${TESTING_OUTPUT_DIR}/perf -p 1 -r 1 -pH 0)
ADD_TEST(NAME PerfusionDerivativesP0R0PH1Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 1 -o ${TESTING_OUTPUT_DIR}/perf -p 0 -r 0 -pH 1)
ADD_TEST(NAME PerfusionDerivativesP1R0PH1Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 1 -o ${TESTING_OUTPUT_DIR}/perf -p 1 -r 0 -pH 1)
ADD_TEST(NAME PerfusionDerivativesP1R0PH0Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 1 -o ${TESTING_OUTPUT_DIR}/perf -p 1 -r 0 -pH 0)
ADD_TEST(NAME PerfusionDerivativesP0R1PH1Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 1 -o ${TESTING_OUTPUT_DIR}/perf -p 0 -r 1 -pH 1)
ADD_TEST(NAME PerfusionDerivativesP0R1PH0Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 1 -o ${TESTING_OUTPUT_DIR}/perf -p 0 -r 1 -pH 0)
ADD_TEST(NAME PerfusionDerivativesEcho5Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 5 -o ${TESTING_OUTPUT_DIR}/perf -p 1 -r 1 -pH 1)
ADD_TEST(NAME PerfusionDerivativesEcho10Test COMMAND PerfusionDerivatives -i ${TESTING_DATA_DIR}/perf/AAAC_PreOp_perf_pp.nii.gz -e 10 -o ${TESTING_OUTPUT_DIR}/perf -p 1 -r 1 -pH 1)

ADD_TEST(NAME PerfusionDerivativesHelpTest COMMAND PerfusionDerivatives -h )
ADD_TEST(NAME PerfusionDerivativesVersionTest COMMAND PerfusionDerivatives -v )

# Perfusion PCA Tests (Incomplete)

ADD_TEST(NAME PerfusionPCAHelpTest COMMAND PerfusionPCA -h )
ADD_TEST(NAME PerfusionPCAVersionTest COMMAND PerfusionPCA -v )

# PopulationAtlases Tests

ADD_TEST(NAME PopulationAtlasesTest COMMAND PopulationAtlases -i ${TESTING_DATA_DIR}/popatlas/atlas/ -l ${TESTING_DATA_DIR}/popatlas/atlaslabelfile.csv -a ${TESTING_DATA_DIR}/popatlas/jakob_stripped_with_cere_lps_256256128.nii.gz -o ${TESTING_OUTPUT_DIR}/popatlas )

ADD_TEST(NAME PopulationAtlasesHelpTest COMMAND PopulationAtlases -h )
ADD_TEST(NAME PopulationAtlasesVersionTest COMMAND PopulationAtlases -v )

# PseudoProgressionEstimator Tests

# Doesn't actually create a new model

# Takes REALLY long
# ADD_TEST(NAME PseudoProgressionEstimatorNewModelTest COMMAND PseudoProgressionEstimator -t 0 -i ${TESTING_DATA_DIR}/pseudo -o ${TESTING_OUTPUT_DIR}/pseudo  )
# ADD_TEST(NAME PseudoProgressionEstimatorExistingModelTest COMMAND PseudoProgressionEstimator -t 1 -i ${TESTING_DATA_DIR}/pseudo -o ${TESTING_OUTPUT_DIR}/pseudo -m ${TESTING_OUTPUT_DIR}/pseudo )

ADD_TEST(NAME PseudoProgressionEstimatorHelpTest COMMAND PseudoProgressionEstimator -h )
ADD_TEST(NAME PseudoProgressionEstimatorVersionTest COMMAND PseudoProgressionEstimator -v )

# Recurrence Estimator Tests

ADD_TEST(NAME RecurrenceEstimatorType0Test COMMAND RecurrenceEstimator -t 0 -i ${TESTING_DATA_DIR}/recurrence/training_samples -o ${TESTING_OUTPUT_DIR}/recurrence/ )
ADD_TEST(NAME RecurrenceEstimatorType1Test COMMAND RecurrenceEstimator -t 1 -i ${TESTING_DATA_DIR}/recurrence/test_samples -o ${TESTING_OUTPUT_DIR}/recurrence/ -m ${TESTING_OUTPUT_DIR}/recurrence/)

ADD_TEST(NAME RecurrenceEstimatorHelpTest COMMAND RecurrenceEstimator -h )
ADD_TEST(NAME RecurrenceEstimatorVersionTest COMMAND RecurrenceEstimator -v )

# SBRT Tests

# Analysis

ADD_TEST(NAME SBRT_AnalysisTest COMMAND SBRT_Analysis -i ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -m ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_testTumor_75.nii.gz -l 1 )
ADD_TEST(NAME SBRT_AnalysisOutputTest COMMAND SBRT_Analysis -i ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -m ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_testTumor_75.nii.gz -l 1 ${TESTING_OUTPUT_DIR}/sbrt/analysis )

ADD_TEST(NAME SBRT_AnalysisHelpTest COMMAND SBRT_Analysis -h )
ADD_TEST(NAME SBRT_AnalysisVersionTest COMMAND SBRT_Analysis -v )

# Lung 
ADD_TEST(NAME SBRT_LungFieldTest COMMAND SBRT_LungField -p ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -c ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -o ${TESTING_OUTPUT_DIR}/sbrt/lung )
ADD_TEST(NAME SBRT_LungFieldMaskTest COMMAND SBRT_LungField -p ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -c ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -o ${TESTING_OUTPUT_DIR}/sbrt/lung -m ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_testTumor_75.nii.gz )

ADD_TEST(NAME SBRT_LungFieldHelpTest COMMAND SBRT_LungField -h )
ADD_TEST(NAME SBRT_LungFieldVersionTest COMMAND SBRT_LungField -v )

# Nodule
ADD_TEST(NAME SBRT_NoduleTest COMMAND SBRT_Nodule -p ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -c ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -m ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_testTumor_75.nii.gz -o ${TESTING_OUTPUT_DIR}/sbrt/nodule )
ADD_TEST(NAME SBRT_NoduleLabelTest COMMAND SBRT_Nodule -p ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -c ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_75.nii.gz -m ${TESTING_DATA_DIR}/sbrt/AAAC0_flair_pp_shrunk_testTumor_75.nii.gz -l 2 -o ${TESTING_OUTPUT_DIR}/sbrt/nodule )

ADD_TEST(NAME SBRT_NoduleHelpTest COMMAND SBRT_Nodule -h )
ADD_TEST(NAME SBRT_NoduleVersionTest COMMAND SBRT_Nodule -v )

# Survival Predictor Tests

ADD_TEST(NAME SurvivalPredictorT0Test COMMAND SurvivalPredictor -t 0 -i ${TESTING_DATA_DIR}/survival -o ${TESTING_OUTPUT_DIR}/survival/ )
ADD_TEST(NAME SurvivalPredictorT1Test COMMAND SurvivalPredictor -t 1 -i ${TESTING_DATA_DIR}/survival -o ${TESTING_OUTPUT_DIR}/survival/ -m ${TESTING_OUTPUT_DIR}/survival/ )

ADD_TEST(NAME SurvivalPredictorHelpTest COMMAND SurvivalPredictor -h )
ADD_TEST(NAME SurvivalPredictorVersionTest COMMAND SurvivalPredictor -v )

# Training Module Tests

ADD_TEST(NAME TrainingModuleTest5Folds COMMAND TrainingModule -n 1 -f ${TESTING_DATA_DIR}/training/features2.csv -l ${TESTING_DATA_DIR}/training/labels2.csv -c 1 -o ${TESTING_OUTPUT_DIR}/training/ -k 5 )
ADD_TEST(NAME TrainingModuleTest10Folds COMMAND TrainingModule -n 1 -f ${TESTING_DATA_DIR}/training/features2.csv -l ${TESTING_DATA_DIR}/training/labels2.csv -c 1 -o ${TESTING_OUTPUT_DIR}/training/ -k 10 )

ADD_TEST(NAME TrainingModuleHelpTest COMMAND TrainingModule -h )
ADD_TEST(NAME TrainingModuleVersionTest COMMAND TrainingModule -v )

# White Stripe Tests 

# White stripe has a bunch of flags that affect a bunch of default values
# since we want to test for "does it work", whatever data we get out we really don't care about and thus
# changing default values theoretically doesn't do anything different to the 
# program execution, so we don't test all of those options.

ADD_TEST(NAME WhiteStripeTest COMMAND WhiteStripe -i ${TESTING_DATA_DIR}/whitestripe/in.nii.gz -o ${TESTING_OUTPUT_DIR}/whitestripe/ )
ADD_TEST(NAME WhiteStripeT1Test COMMAND WhiteStripe -t1=0 -i ${TESTING_DATA_DIR}/whitestripe/in.nii.gz -o ${TESTING_OUTPUT_DIR}/whitestripe/ )

ADD_TEST(NAME WhiteStripeHelpTest COMMAND WhiteStripe -h )
ADD_TEST(NAME WhiteStripeVersionTest COMMAND WhiteStripe -v )

##

# Utilities Test-1
# CaPTk Utilities not created
# ADD_TEST(NAME UtilitiesTestHistoMatch COMMAND CaPTk_Utilities -i ${TESTING_DATA_DIR}/AAAC0_flair_pp_shrunk.nii.gz -o ${TESTING_OUTPUT_DIR}/histoMatch.nii.gz -hm 40 -hl 120 )

# # buffer test
# ADD_TEST(NAME BufferTest COMMAND ${TEST_EXE_NAME} --bufferTest "random" )

# # geodesic test
# ADD_TEST(NAME GeodesicTest COMMAND ${TEST_EXE_NAME} --geodesic "${TESTING_DATA_DIR}/Geodesic/" )

# # EGFRvIII test
# ADD_TEST(NAME EGFRvIIITest COMMAND ${TEST_EXE_NAME} --egfrviii "${TESTING_DATA_DIR}/EGFRvIII/" )
