/**
\file  fMainWindow.h

\brief Declaration of fMainWindow class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _fMainWindow_h_
#define _fMainWindow_h_

#include "NiftiDataManager.h"
#include "RecurrenceEstimator.h"
#include "PseudoProgressionEstimator.h"
#include "SurvivalPredictor.h"
#include "EGFRvIIIIndexPredictor.h"
#include "ui_fMainWindow.h"
#include "OutputWritingManager.h"
#include "PreprocessingPipelineClass.h"
#include "FeatureExtractionClass.h"
#include "fHelpTutorial.h"
#include "PopulationAtlases.h"
#include "ImagingSubtypePredictor.h"
#include "MolecularSubtypePredictor.h"
#include "FetalBrain.h"
#include "SlicerManager.h"
#include "itkVTKImageImport.h"
#include "vtkImageExport.h"
#include "vtkImageData.h"

#include "fTumorPanel.h"
#include "fImagesPanel.h"
#include "fDrawingPanel.h"
#include "fTrainingDialog.h"
#include "fFeaturePanel.h"
#include "fRecurrenceDialog.h"
#include "fPseudoProgressionDialog.h"
#include "fRegistrationDialog.h"
#include "fPreprocessingDialog.h"
#include "fSurvivalDialog.h"
#include "fEGFRvIIIDialog.h"
#include "fSkullStripDialog.h"
#include "fPerfusionMeasuresDialog.h"
#include "fPerfusionAlignmentDialog.h"
#include "fDiffusionMeasuresDialog.h"
#include "fPCADialog.h"
#include "fHistoMatchDialog.h"
#include "fWhiteStripeDialog.h"
#include "fDirectionalityDialog.h"
#include "fPopulationAtlasDialog.h"
#include "fImagingSubtypeDialog.h"
#include "fMolecularSubtypeDialog.h"
#include "fDCM2NIfTI.h"
#include "fDeepMedicDialog.h"
#include "fTexturePipelineDialog.h"
#include "fDeepMedicNormDialog.h"
#include "fFetalBrain.h"
#include "fSBRTNoduleDialog.h"
#include "fSBRTAnalysisDialog.h"
#include "fBiasCorrectionDialog.h"
#include "fBraTSSegmentation.h"

#include <atomic>

#include "GeodesicTrainingCaPTkApp.h"

#include <QMessageBox>

#include "itkJoinSeriesImageFilter.h"
#include "itkExtractImageFilter.h"

#include "QVTKOpenGLWidget.h"
#include <QScopedPointer>
#include "vtkGenericOpenGLRenderWindow.h"
#include "fBottomImageInfoTip.h"

#include "ApplicationDownloadManager.h"
#include "yaml-cpp/node/node.h"

class SlicerManager;
class Slicer;
class SimpleImageManager;
class fHelpDialog;
class PreferencesDialog;
class SystemInformationDisplayWidget;

#define USE_PROCESSDIALOG

enum TAB_TYPE
{
  TAB_IMAGES, TAB_TUMOR, TAB_DRAW, TAB_FEWATURE
};

typedef itk::Image< short, 3 > GenericImage;

// multiLabel
enum DRAW_MODE
{
  DRAW_MODE_LABEL_1 = 1, DRAW_MODE_LABEL_2, DRAW_MODE_LABEL_3, DRAW_MODE_LABEL_4, DRAW_MODE_LABEL_5, DRAW_MODE_LABEL_6,
  DRAW_MODE_LABEL_7, DRAW_MODE_LABEL_8, DRAW_MODE_LABEL_9
};

enum SHAPE_MODE
{
  SHAPE_MODE_NONE = 0, SHAPE_MODE_ERASER, SHAPE_MODE_FREE_HAND, SHAPE_MODE_LINE, SHAPE_MODE_RECTANGLE, SHAPE_MODE_CIRCLE, SHAPE_MODE_FILL, SHAPE_MODE_SPHERE
};

//TBD move all this to util class 
template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines(ITK_Exporter &exporter, VTK_Importer* importer)
{
  importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  importer->SetSpacingCallback(exporter->GetSpacingCallback());
  importer->SetOriginCallback(exporter->GetOriginCallback());
  importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  importer->SetCallbackUserData(exporter->GetCallbackUserData());
}
template <typename VTK_Exporter, typename ITK_Importer>
void ConnectPipelines(VTK_Exporter* exporter, ITK_Importer &importer)
{
  importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  importer->SetSpacingCallback(exporter->GetSpacingCallback());
  importer->SetOriginCallback(exporter->GetOriginCallback());
  importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  importer->SetCallbackUserData(exporter->GetCallbackUserData());
}

//! Converts a VTK-Image to ITK-Image on the basis of the pixel type definition of the former and latter
template<class InputPixelType, class OutputPixelType, int VImageDimension>
typename itk::Image<OutputPixelType, VImageDimension>::Pointer convertVtkToItk(vtkSmartPointer<vtkImageData> vtkInput)
{
  typedef itk::Image<InputPixelType, VImageDimension> InputImageType;
  typedef itk::Image<OutputPixelType, VImageDimension> OutputImageType;
  typedef itk::VTKImageImport<InputImageType> InputImageImportType;

  vtkSmartPointer<vtkImageExport> inputImageExporter = vtkImageExport::New();
#if VTK_MAJOR_VERSION <= 5
  inputImageExporter->SetInput(vtkInput);
#else
  inputImageExporter->SetInputData(vtkInput);
#endif
  typename InputImageImportType::Pointer inputImageImporter = InputImageImportType::New();

  ConnectPipelines(inputImageExporter.GetPointer(), inputImageImporter);

  typename InputImageType::Pointer inputImage = const_cast<InputImageType*>(inputImageImporter->GetOutput());
  inputImage->Update();
  if (inputImage.IsNull()) {
    std::cerr << "vtkInput is invalid." << std::endl;
  }

  typedef itk::CastImageFilter< InputImageType, OutputImageType > CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(inputImage);
  castFilter->Update();

  typename OutputImageType::Pointer outputImage = castFilter->GetOutput();
  outputImage->Update();
  if (outputImage.IsNull()) {
    std::cerr << "inputImage is invalid." << std::endl;
    return NULL;
  }

  return outputImage;
}

//! Converts a VTK-Image to ITK-Image on the basis of the pixel type definition of the latter
template<class OutputPixelType, int VImageDimension>
typename itk::Image<OutputPixelType, VImageDimension>::Pointer convertVtkToItk(vtkSmartPointer<vtkImageData> vtkInput) //TBD move to cbica util 
{
  std::string InputPixelType = vtkInput->GetScalarTypeAsString();
  if (InputPixelType == "short")
  {
    return convertVtkToItk<short, OutputPixelType, VImageDimension>(vtkInput);
  }
  else if (InputPixelType == "unsigned short")
  {
    return convertVtkToItk<unsigned short, OutputPixelType, VImageDimension>(vtkInput);
  }
  else if (InputPixelType == "char")
  {
    return convertVtkToItk<char, OutputPixelType, VImageDimension>(vtkInput);
  }
  else if (InputPixelType == "unsigned char")
  {
    return convertVtkToItk<unsigned char, OutputPixelType, VImageDimension>(vtkInput);
  }
  else if (InputPixelType == "int")
  {
    return convertVtkToItk<int, OutputPixelType, VImageDimension>(vtkInput);
  }
  else if (InputPixelType == "unsigned int")
  {
    return convertVtkToItk<unsigned int, OutputPixelType, VImageDimension>(vtkInput);
  }
  else if (InputPixelType == "double")
  {
    return convertVtkToItk<double, OutputPixelType, VImageDimension>(vtkInput);
  }
  else if (InputPixelType == "float")
  {
    return convertVtkToItk<float, OutputPixelType, VImageDimension>(vtkInput);
  }
  else {
    std::cerr << "Error, input pixel type : " << InputPixelType << " unknown !" << std::endl;
  }
  return NULL;
}

/**
\class fMainWindow

\brief This is the main UI class for CaPTk
*/
#if __GNUC__
#pragma GCC visibility push(hidden)
#endif
class fMainWindow : public QMainWindow, private Ui::fMainWindow
{
  Q_OBJECT

private:
  fRecurrenceDialog			recurrencePanel;
  fPseudoProgressionDialog pseudoPanel;
  fPopulationAtlasDialog	atlasPanel;
  fRegistrationDialog		registrationPanel;
  fPreprocessingDialog	preprocessingPanel;
  fSurvivalPredictor survivalPanel;
  fEGFRvIIIPredictor egfrv3Panel;
  fMolecularSubtypePredictor msubtypePanel;
  fImagingSubtypePredictor isubtypePanel;
  fFetalBrain fetalbrainpanel;
  fSBRTNoduleDialog nodulePanel;
  fSBRTAnalysisDialog analysisPanel;
  fBiasCorrectionDialog biascorrectionPanel;

  fSkullStripper skullStrippingPanel;
  fPCADialog pcaPanel;
  fTrainingSimulator trainingPanel;
  fPerfusionEstimator perfmeasuresPanel;
  fPerfusionAligner perfalignPanel;
  fDiffusionEstimator diffmeasuresPanel;
  fDCM2NIfTIConverter dcmConverter;
  fDeepMedicDialog deepMedicDialog;
  fTexturePipelineDialog texturePipelineDialog;
  fHistoMatcher histoMatchPanel;
  fDeepMedicNormalizer deepMedicNormPanel;
  fWhiteStripeObj whiteStripeNormalizer;
  fBraTSSegmentation bratsPipelineDialog;
  fDirectionalityDialog directionalityEstimator;
  PreferencesDialog *preferenceDialog;
  ApplicationDownloadManager appDownloadMngr;
  SystemInformationDisplayWidget *sysinfowidget;
  

  fDrawingPanel *drawingPanel;
  fFeaturePanel *featurePanel;
  fImagesPanel *imagesPanel;
  fBottomImageInfoTip *infoPanel;
  fTumorPanel *tumorPanel;

  //-------------menu-----------
  QMenuBar *menubar;
  QMenu* menuFile;
  QMenu* menuLoadFile;
  QMenu* menuSaveFile;
  QMenu* menuExit;
  QMenu* menuLoadFileDicom;
  QMenu* menuLoadFileNifti;
  QMenu* menuDownload;
  QMenu* menuApp;
  QMenu* menuPreprocessing;
  QMenu* menuDeepLearning;
  QMenu* menuHelp;

  QAction *help_discussion;
  QAction *helpMenu_download;
  QAction *help_forum;
  QAction *help_bugs;
  QAction *help_features;
  QAction *help_systeminformation;

  //-------------actions-------------

  QAction *actionLoad_Recurrence_Images;
  QAction *actionLoad_Nifti_Images;
  QAction *actionLoad_Nifti_ROI;
  QAction *actionLoad_Dicom_Images;
  QAction *actionPreferences;


  QAction *actionSave_Nifti_Images;
  QAction *actionSave_Dicom_Images;
  QAction *actionSave_ROI_Images;
  QAction *actionSave_ROI_Dicom_Images;

  QAction *actionHelp_Interactions;
  QAction *actionSave_Images;
  QAction *actionAbout;
  QAction *actionExit;
  QAction *actionModelLibrary;

  QAction *actionAppEGFR;
  QAction *actionAppRecurrence;
  QAction *actionAppGeodesic;
  QAction *actionAppGeodesicTraining;

  // obtain list from CMake variables using populateStringListInMenu() function
  std::vector< std::string >
    m_nativeApps, // native CPP applications
    m_preprocessApps, // native pre-processing routines
    m_pyCLIApps, // python command line applications
    m_pyGUIApps; // python graphical applications

  std::map< std::string, std::string > m_allNonNativeApps;

  //! QVTK OpenGL Widgets
  QScopedPointer<QVTKOpenGLWidget> SaggitalViewWidget;
  QScopedPointer<QVTKOpenGLWidget> AxialViewWidget;
  QScopedPointer<QVTKOpenGLWidget> CoronalViewWidget;

  //! renderwindows
  vtkSmartPointer< vtkGenericOpenGLRenderWindow> SaggitalRenWin;
  vtkSmartPointer< vtkGenericOpenGLRenderWindow> AxialRenWin;
  vtkSmartPointer< vtkGenericOpenGLRenderWindow> CoronalRenWin;

  QHBoxLayout* bottomLayout;
    
  YAML::Node m_downloadLinks; //! structure to save download links

  /**
  \struct ActionAndName

  \brief This is a helper struct to tie an action with its name as a std::string
  */
  struct ActionAndName
  {
    QAction* action;
    std::string name;
  };

  //! Wrap to ensure previous functionality doesn't break
  std::vector< ActionAndName > populateStringListInMenu(const std::string &inputList, QMainWindow* inputFMainWindow, QMenu* menuToPopulate, std::string menuAppSubGroup, bool ExcludeGeodesic);

  /**
  \brief Takes a list of application variables from CMake defines and put it in specified window and menu

  \param inputList The list obtained from CMake variable which is added to cache
  \param inputFMainWindow The current fMainWindow from which the QActions need to inherit
  \param menuToPopulate The QMenu in which the QActions need to be populated *visualizationInputImagesLabel
  \return A vector of ActionAndName structs which ties a QAction to the corresponding name from inputList
  **/
  std::vector< ActionAndName > populateStringListInMenu(const std::vector< std::string > &vectorOfInputs, QMainWindow* inputFMainWindow, QMenu* menuToPopulate, std::string menuAppSubGroup, bool ExcludeGeodesic);

  //! check if input files also include directories
  bool hasDirectories(QStringList &lst, int &nDirs);

  // initialize vectors of Actions and Names so that the process can be automated and the QAction is tied to its corresponding Name
  std::vector< ActionAndName >
    vectorOfGBMApps, // GBM-specific applications
    vectorOfBreastApps, // breast-specific applications
    vectorOfLungApps, // lung-specific applications
    vectorOfSegmentationApps, // the segmentation applications
    vectorOfMiscApps, // the rest
    vectorOfPreprocessingActionsAndNames, // for preprocessing algorithms
    vectorOfDeepLearningActionsAndNames; // for deep learning applications

  QTableWidget * m_imagesTable;
  QTableWidget * m_nonVisImagesTable;

public:
  //! Default constructor
  fMainWindow();

  //! Default destructor
  ~fMainWindow();

  //! Check if a valid mask is defined for the current instance of the 
  bool isMaskDefined();

  //! Set Comparison Mode
  void SetComparisonMode(bool mode);

  //! Get Comparison Mode status
  bool GetComparisonMode();

  /*
  \brief Change direction cosine of image to identity
  \param image The image whose direction is altered
  */
  template< class TImageType = ImageTypeFloat3D >
  typename TImageType::Pointer ChangeImageDirectionToIdentity(const typename TImageType::Pointer image)
  {
    auto infoChangeFilter = itk::ChangeInformationImageFilter< TImageType >::New();
    typename TImageType::DirectionType direction_new;
    direction_new.SetIdentity();
    infoChangeFilter->SetInput(image);
    infoChangeFilter->SetOutputDirection(direction_new);
    infoChangeFilter->ChangeDirectionOn();
    infoChangeFilter->Update();
    return infoChangeFilter->GetOutput();
  }

  /**
  \brief Load images into memory

  \param filenames Vector of image filenames which need to be loaded (comes from "Load" dialog)
  \param imagetype_int Either NIfTI or DICOM
  \param bSkipDup Skip duplicates, defaults to true
  */
  void LoadSlicerImages(const std::string &fileName, const int &imagetype_int, bool bSkipDup = true);

  /**
  \brief Load non-viewing images into memory

  \param directoryname Directory name in which the non-viewing images are present
  \param imagetype_int Either NIfTI or DICOM
  \param imagesubtype Modality of the image (T1, T2, ...)
  */
  void LoadNonViewingImages(const std::string &directoryname, const int &imagetype_int, const int &imagesubtype);

  /**
  \brief Initialize drawing mask on the image
  */
  void InitMask(vtkImageData* image);
   
  /**
  \brief Initializer for the entire rendering pipeline
  */
  void InitDisplay();

  /**
  \brief Initializes difference slicer views of a single image
  */
  void InitSlicers();

  /**
  \brief Render the sliders to adjust the view of the image slices.

  Initializes it at the middle of the whole image.

  \param slicer If many images are present, this keeps a count of which image is being rendered
  \param window The window number (goes from 1 to 3)
  */
  void DisplaySliders(int slicer, int window);

  /**
  \brief Get image index from the visualization table
  */
  int GetSlicerIndexFromItem(QTableWidgetItem* item);

  /**
  \brief Get the corresponding table item from the image
  */
  QTableWidgetItem* GetItemFromSlicerManager(SlicerManager* sm);

  /**
  \brief Construct near and far indeces from the initialized mask
  */
  template <class TImageType>
  typename TImageType::Pointer GetImageWithLabels(std::vector<double> labels, typename TImageType::Pointer inputimage)
  {
    auto img = convertVtkToItk<float, 3>(mSlicerManagers[0]->mMask);

    typename TImageType::Pointer output = TImageType::New();
    output->CopyInformation(img);
    output->SetRequestedRegion(img->GetLargestPossibleRegion());
    output->SetBufferedRegion(img->GetBufferedRegion());
    output->Allocate();
    output->FillBuffer(0);

    typedef itk::ImageRegionIteratorWithIndex< TImageType > IteratorType;
    IteratorType maskIt(img, img->GetLargestPossibleRegion());
    IteratorType outputIt(output, output->GetLargestPossibleRegion());
    maskIt.GoToBegin();
    outputIt.GoToBegin();
    while (!maskIt.IsAtEnd())
    {
      for (int j = 0; j < labels.size(); j++)
      {
        if (maskIt.Get() == labels[j])
        {
          outputIt.Set(CAPTK::VOXEL_STATUS::ON);
        }
      }
      ++maskIt;
      ++outputIt;
    }
    return output;
  }

  /**
  \brief Rescales image intensity in the range of 0-255
  */
  ImageTypeFloat3D::Pointer RescaleImageIntensity(ImageTypeFloat3D::Pointer image);


  /*
  \brief This function is used to load parameters from the command line

  \param files The input image file(s)
  \param comparisonMode Flag that enables/disables comparison mode
  \param maskImage The mask image corresponding to the input image(s)
  \param maskOpacity The opacity of the loaded ROI labels
  \param tumorPointFile The tumor points file (containing radius information)
  \param tissuePointFile The tissue points file (containing no radius information)
  */
  void loadFromCommandLine(std::vector< QString > files, bool comparisonMode, const std::string &maskImage = "", const float maskOpacity = 1.0,
    const std::string &tumorPointFile = "", const std::string &tissuePointFile = "", bool firstRun = false);

signals:
  void SelectedImageHasChanged(SlicerManager *);
  void LandmarksFocused(bool bFocused);
  void SeedPointsFocused(bool bFocused);
  void TissuePointsFocused(bool bFocused);

public slots:

	//! apply mask to loaded images
	void OnApplyMask();

	//!display Preferences dialog
	void OnPreferencesMenuClicked();

  //! set Z slice position on image info panel
  void SetImageInfoZSlicePosition(int zslice);

  //! set voxel intensity value at cursor position on image info panel
  void SetImageInfoIntensityValue(double value);

  //! slot on movement of slider in comparison mode
  void OnSliderMovedInComparisonMode(int);

  /**
  \brief Drag event initialization. Can accept events emitted from widgets.
  */
  void dragEnterEvent(QDragEnterEvent *event);

  /**
  \brief Drop Event parsing slot. Can accept events emitted from widgets.
  */
  void dropEvent(QDropEvent *event);

  /**
  \brief Updates draw mode when drawing panel changes
  */
  void updateDrawMode(int mode = -1);

  /**
  \brief Updates number of near and far points in the table
  */
  void UpdateNumberOfPointsInTable();

  /**
  \brief Initializes preset combo-box with the default options
  */
  void SetPresetComboBox();

  /**
  \brief Main function that performs EGFRvIII estimation on the displayed image
  */
  void StartEGFREstimate();

  /**
  \brief Get the mask as an ITK image as it is currently displayed
  */
  ImageTypeFloat3D::Pointer getMaskImage();
  
  /**
  \brief Function that performs segmentation and UI measurements of Fetal Ventriculomegaly application
  */
  void FetalBrain_SkullStripfunc();
  
  /**
  \brief Function that performs prediction of Fetal Ventriculomegaly application
  */
  void FetalBrain_Predict();
  
  /**
  \brief Function that performs new model training for Fetal Ventriculomegaly application
  */
  void FetalBrain_TrainNewModel(const std::string &directory, const std::string &outputdirectory);

  /**
  \brief get images loaded and their file names
  */
  std::vector<ImageTypeFloat3D::Pointer> getLodedImages(std::vector<std::string> &fileNames, std::vector<std::string> &modality, bool onlySelected = false);

  /**
  \brief Main function that estimates recurrence on the displayed test subject
  param outputdirectory The directory where recurrence map will be written
  param cbT1Data Whether T1 data is present or not
  param cbT1ceData Whether T1CE data is present or not
  param cbT2Data Whether T2 data is present or not
  param cbT2FlairData Whether T2-Flair data is present or not
  param cbDTIData Whether DTI data is present or not
  param cbPerfData Whether Perfusion data is present or not
  param cbDistData Whether Distance feature need to be used or not
  */
  void StartRecurrenceEstimate(const std::string &outputdirectory, bool cbT1Data, bool cbDTIData, bool cbPerfData, bool cbDistData);

  /**
  \brief Load the subjects for an existing Reccurence model estimate
  param outputdirectory The directory where recurrence map will be written
  param outputdirectory The output directory
  param modeldirectory The input model directory
  param cbT1Data Whether T1 data is present or not
  param cbPerfData Whether Perfusion data is present or not
  param cbT2Data Whether T2 data is present or not
  param cbDistData Whether Distance data is present or not
  param cbDTIData Whether DTI data is present or not
  */
  void LoadedSubjectExistingRecurrenceEstimate(const std::string &outputdirectory, const std::string &modeldirectory, bool cbT1Data, bool cbDTIData, bool cbPerfData, bool cbDistData);

  /**
  \brief Main function that estimates recurrence on test subjects by using an existing model
  param modeldirectory The directory where model related files are stored
  param inputdirectory The directory where test data is stored
  param outputdirectory The directory where output data will be stored

  param cbT1Data Whether T1 data is present or not
  param cbT1ceData Whether T1CE data is present or not
  param cbT2Data Whether T2 data is present or not
  param cbT2FlairData Whether T2-Flair data is present or not
  param cbDTIData Whether DTI data is present or not
  param cbPerfData Whether Perfusion data is present or not
  param cbDistData Whether Distance feature need to be used or not
  */
  void RecurrenceEstimateOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::string &outputdirectory, bool cbConventionalData, bool cbDTIData, bool cbPerfData, bool cbDistData);

  /**
  \brief Main function that trains a model on many training subjects
  param directory The directory that contains training data
  param outputdirectory The directory where trained model related data will be stored
  param cbT1Data Whether T1 data is present or not
  param cbT1ceData Whether T1CE data is present or not
  param cbT2Data Whether T2 data is present or not
  param cbT2FlairData Whether T2-Flair data is present or not
  param cbDTIData Whether DTI data is present or not
  param cbPerfData Whether Perfusion data is present or not
  param cbDistData Whether Distance feature need to be used or not
  */
  void TrainNewModelOnGivenData(const std::string &directory, const std::string &outputdirectory, bool cbConvData, bool cbDTIData, bool cbPerfData, bool cbDistData);

  /**
  \brief Main function that estimates pseudoprogression on the displayed test subject
  param outputdirectory The directory where recurrence map will be written
  param cbT1Data Whether T1 data is present or not
  param cbT1ceData Whether T1CE data is present or not
  param cbT2Data Whether T2 data is present or not
  param cbT2FlairData Whether T2-Flair data is present or not
  param cbDTIData Whether DTI data is present or not
  param cbPerfData Whether Perfusion data is present or not
  param cbDistData Whether Distance feature need to be used or not
  */
  void LoadedSubjectExistingPseudoprogressionEstimate(const std::string &outputdirectory, const std::string &modeldirectory, bool cbT1Data, bool cbDTIData, bool cbPerfData, bool cbDistData);

  /**
  \brief Main function that estimates PseudoProgression on test subjects by using an existing model
  param modeldirectory The directory where model related files are stored
  param inputdirectory The directory where test data is stored
  param outputdirectory The directory where output data will be stored

  param cbT1Data Whether T1 data is present or not
  param cbT1ceData Whether T1CE data is present or not
  param cbT2Data Whether T2 data is present or not
  param cbT2FlairData Whether T2-Flair data is present or not
  param cbDTIData Whether DTI data is present or not
  param cbPerfData Whether Perfusion data is present or not
  param cbDistData Whether Distance feature need to be used or not
  */
  void PseudoprogressionEstimateOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::string &outputdirectory, bool cbConventionalData, bool cbDTIData, bool cbPerfData, bool cbDistData);

  /**
  \brief Main function that trains a model on many training subjects
  param directory The directory that contains training data
  param outputdirectory The directory where trained model related data will be stored
  param cbT1Data Whether T1 data is present or not
  param cbT1ceData Whether T1CE data is present or not
  param cbT2Data Whether T2 data is present or not
  param cbT2FlairData Whether T2-Flair data is present or not
  param cbDTIData Whether DTI data is present or not
  param cbPerfData Whether Perfusion data is present or not
  param cbDistData Whether Distance feature need to be used or not
  */
  void TrainNewPseudoprogressionModelOnGivenData(const std::string &directory, const std::string &outputdirectory, bool cbConvData, bool cbDTIData, bool cbPerfData, bool cbDistData);

  /**
  \brief Survival analysis using Existing model

  \param modeldirectory The directory where the existing model files are stored
  \param inputdirectory The input subject directory (pre-sorted)
  \param outputdirectory The output directory to save the data
  */
  void CallForSurvivalPredictionOnExistingModelFromMain(const std::string modeldirectory, const std::string inputdirectory, const std::string outputdirectory);

  /**
  \brief Create new model for Survival analysis

  \param inputdirectory The input subjects directory (pre-sorted)
  \param outputdirectory The output directory to save the data
  */
  void CallForNewSurvivalPredictionModelFromMain(const std::string inputdirectory, const std::string outputdirectory);
  
  /**
  \brief EGFRvIII estimation using Existing model

  \param modeldirectory The directory where the existing model files are stored
  \param inputdirectory The input subject directory (pre-sorted)
  \param outputdirectory The output directory to save the data
  */
  void CallForEGFRvIIIPredictionOnExistingModelFromMain(const std::string modeldirectory, const std::string inputdirectory, const std::string outputdirectory);

  /**
  \brief Create new model for EGFRvIII prediction

  \param inputdirectory The input subjects directory (pre-sorted)
  \param outputdirectory The output directory to save the data
  */
  void CallForNewEGFRvIIIPredictionModelFromMain(const std::string inputdirectory, const std::string outputdirectory);
  
  /**
  \brief MolecularSubtype using Existing model

  \param modeldirectory The directory where the existing model files are stored
  \param inputdirectory The input subject directory (pre-sorted)
  \param outputdirectory The output directory to save the data
  */
  void CallForMolecularSubtypePredictionOnExistingModelFromMain(const std::string modeldirectory, const std::string inputdirectory, const std::string outputdirectory);

  /**
  \brief Create new model for MolecularSubtype

  \param inputdirectory The input subjects directory (pre-sorted)
  \param outputdirectory The output directory to save the data
  */
  void CallForNewMolecularSubtypePredictionModelFromMain(const std::string inputdirectory, const std::string outputdirectory);

  /**
  \brief Call the Skull Stripping module from ITK with the inputs

  \param referenceAtlas The reference atlas (defaults to SRI24 atlas)
  \param referenceMask The reference mask (defaults to SRI24 mask)
  \param inputImageFile The input image file
  \param outputImageFile The output file to save
  */
  void CallImageSkullStripping(const std::string referenceAtlas, const std::string referenceMask, const std::string inputImageFile, const std::string outputImageFile);

  /**
  \brief Call the Perfusion alignment application with the inputs

  \param echotime The echo time of the input Perfusion image
  \param before Number of time-points before the standard perfusion curve that we want our image to be aligned with.
  \param after Number of time-points after the standard perfusion curve that we want our image to be aligned with.
  \param inputfilename The input perfusion image file name
  \param inputt1cefilename The input T1-Gd file name
  \param inputdicomfilename The input DICOM slide
  \param outputFolder The output folder to write all results
  */
  void CallPerfusionAlignmentCalculation(const double echotime, const int before, const int after, const std::string inputfilename, const std::string inputt1cefilename, std::string outputFolder);

  /**
  \brief Call the Perfusion Measures application with the inputs

  \param TE The echo time of the input Perfusion image
  \param rcbv Flag that enables RCBV calculation
  \param psr Flag that enables PSR calculation
  \param ph Flag that enables PH calculation
  \param inputfile The input DSC-MRI image
  \param outputFolder The output folder to write all results
  */
  void CallPerfusionMeasuresCalculation(const bool rcbv, const bool psr, const bool ph, const std::string inputfile, const std::string outputFolder);

  /**
  \brief Call the Diffusion Measures application with the inputs

  \param inputImage The input DSC-MRI image
  \param maskImage The mask file
  \param BValFile The BVal file
  \param BVecFile The BVec file
  \param ax Flag that enables AX calculation
  \param FA Flag that enables FA calculation
  \param RAD Flag that enables RAD calculation
  \param TR Flag that enables TR calculation
  \param outputFolder The output folder to write all results
  */
  void CallDiffusionMeasuresCalculation(const std::string inputImage, const std::string maskImage, const std::string BValFile, const std::string BVecFile, const bool ax, const bool fa, const bool rad, const bool tr, const std::string outputFolder);

  /**
  \brief Call the Diffusion Measures application with the inputs

  \param featuresfile The input features 
  \param targetfile The labels, rows should be same as those in featuresfile
  \param outputFolder The output folder to write all results
  \brief modeldirectory The directory with saved model
  \param classifier The classifer type
  \param conf The configuration type
  \param folds The number of folds
  */
  void CallTrainingSimulation(const std::string featuresfile, const std::string targetfile, const std::string outputFolder, const std::string modeldirectory, int classifier, int conf, int folds);

  /**
  \brief Call the PCA calculation application with the inputs

  \param modeldirectory The trained model directory
  \param outputdirectory The folder to save results
  */
  void PCAEstimateOnExistingModel(const std::string &modeldirectory, const std::string &inputdirectory, const std::string &outputdirectory);
  
  /**
  \brief Call the PCA calculation application with the inputs

  \param inputdirectory The input folder containing the image(s)
  \param outputdirectory The folder to save the trained model
  */
  void TrainNewPCAModelOnGivenData(const std::string &inputdirectory, const std::string &outputdirectory);
  /**
  \brief Call DCM2NII for DICOM conversion and load the image into CaPTk

  \param firstImageInSeries First image in the DICOM series (assumes the folder has a single DICOM series)
  \param outputFileName The output NIfTI file
  */
  void CallDCM2NIfTIConversion(const std::string firstImageInSeries, bool loadAsImage);

  /**
  \brief Call DCM2NII for DICOM conversion

  \param firstImageInSeries First image in the DICOM series (assumes the folder has a single DICOM series)
  \param outputFileName The output NIfTI file
  */
  void CallDCM2NIfTIConversion(const std::string firstImageInSeries, const std::string outputFileName);

  /**
  \brief Call the Deep Medic Segmentation dialog
  */
  void CallDeepMedicSegmentation(const std::string modelDirectory, const std::string outputDirectory);

  /**
  \brief Call the breast texture pipeline
  */
  void CallTexturePipeline(const std::string outputDirectory);

  /**
  \brief Call Histogram Matching module of ITK
  */
  void CallImageHistogramMatching(const std::string referenceImage, const std::string inputImageFile, const std::string outputImageFile);

  /**
  \brief Call BraTS Pipeline application
  */
  void CallBraTSPipeline(const std::string t1ceImage, const std::string t1Image, const std::string t2Image, const std::string flImage, const std::string outputDir);

  /**
  \brief Call Histogram Matching module of ITK
  */
  void CallLabelValuesChange(const std::string oldValues, const std::string newValues);

  /**
  \brief Call Deep Medic Normalization module
  */
  void CallImageDeepMedicNormalizer(const std::string inputImage, const std::string maskImage, const std::string outputImageFile,
    const std::string quantLower, const std::string quantUpper,
    const std::string cutoffLower, const std::string cutoffUpper, bool wholeImageMeanThreshold);

  //! Call the Directionality Estimator
  void CallDirectionalityEstimator(const std::string roi1File, const std::string roi2File, const std::string outputDir);

  void CallWhiteStripe(double twsWidth, int sliceStartZ, int sliceStopZ, int tissuesMax, double smoothMax, double smoothDelta, int histSize,
    bool T1Image, const std::string outputFileName);

  /**
  \brief Generate population atlas
  */
  void CallGeneratePopualtionAtlas(const std::string inputdirectory, const std::string inputatlas, const std::string outputImageFile);

  /**
  \brief Generete SBRT Nodule
  */
  void CallSBRTNodule(const std::string seedImage, const int labelValue);

  /** 
  \brief Generate and load bias-corrected image
  param correctionType string to determine bias correction method. Accepts n3 or n4 (case-insensitive)
  param saveFileName where to save the output
  param bias_splineOrder
  param bias_otsuBins
  param bias_maxIterations
  param bias_fittingLevels 
  param bias_filterNoise 
  param bias_fwhm full width at half-maximum 
  */
  void CallBiasCorrection(const std::string correctionType, QString saveFileName,
      int bias_splineOrder, int bias_otsuBins, int bias_maxIterations, int bias_fittingLevels,
      float bias_filterNoise, float bias_fwhm);

  /**
  \brief Function that updates the co-ordinates (X and Y) of border
  param startX starting X co-ordinate
  param startY starting Y co-ordinate
  param endX ending X co-ordinate
  param endY ending Y co-ordinate
  */

  void UpdateBorderWidget(double startX, double startY, double endX, double endY);

  /**
  \brief Function that updates the co-ordordinatesicates (Z) of border
  param startX starting Z co-ordinate
  param startY starting Z co-ordinate
  */
  void UpdateBorderWidget(double startZ, double endZ);

  /**
  \brief Adds actions in the action log
  param points The voxels corresponding to action
  */
  void UpdateActionQ(const QVariantList& list)
  {
    std::vector< PointVal > points;
    for (int i = 0; i < list.size(); i++)
    {
      points.push_back(qvariant_cast<PointVal>(list[i]));
    }
    UpdateAction(points);
  }
  void UpdateAction(std::vector< PointVal > points = std::vector< PointVal >());

  /**
  \brief Brief description about the software
  */
  void about();

  /**
  \brief Help for interactions
  */
  void help_Interactions();

  //! Don't show tutorial again
  void skipTutorial(bool flag)
  {
    m_skipTutorialOnNextRun = flag;
  }
 
  /**
  \brief Help for downloading Sample Data
  */
  void help_Download(QAction* action);

  /**
  \brief open model library webpage
  */
  void OpenModelLibrary();

  /**
 \brief system information menu click 
 */
  void OnSystemInformationMenuClicked();

  /**
  \brief Get contextual help 

  \param startPage The starting page for the web engine view
  */
  void help_contextual(const std::string startPage);

  /**
  \brief Open Nifti image functionality. Shows dialog to select nifti files
  */
  void openImages(QStringList files = QStringList(), bool callingFromCmd = false);

  /**
  \brief Open Dicom image functionality. Shows dialog to select Dicom files
  */
  void openDicomImages(QString dir);

  /**
  \brief Function called when the sliders of axial view is changed
  */
  void AxialViewSliderChanged();

  /**
  \brief Function called when the sliders of coronal view is changed
  */
  void CoronalViewSliderChanged();

  /**
  \brief Function called when the sliders of saggital view is changed
  */
  void SaggitalViewSliderChanged();

  /**
  \brief Closing viewing image by pressing X in front of the image
  param item The current selected item of the table
  */
  void CloseImage(QTableWidgetItem* item);

  /**
  \brief Closing non-viewing image by pressing X in front of the image
  param item The current selected item of the table
  */
  void CloseNonViewingDTIImage(QTableWidgetItem* item);

  //! Clear the mask image
  void clearMask(int label = -1);

  //! Make a drawing "stroke" - which can be fed into the "undo" system
  void makeStroke(std::vector<itk::Image<short, 3>::IndexType>& indices, const int value);

  //! Get selected drawing label
  int getSelectedDrawLabel()
  {
    return drawingPanel->getSelectedDrawLabel();
  }
  
  //! Get selected drawing size
  int getSelectedDrawSize()
  {
    return drawingPanel->getSelectedDrawSize();
  }

  /**
  \brief Save near/far drawing in Nifti format
  */
  void SaveDrawing();

  /**
  \brief Save near/far drawing in DICOM format
  */
  void SaveDicomDrawing();

  /**
  \brief Save initial seed drawing in Nifti format
  */
  void SaveSeedDrawing();

  /**
  \brief Load annotated ROI from qt file accept box
  */
  void LoadDrawing();

  /**
  \brief Load annotated ROI from filename
  */
  void LoadDrawing(const std::string &maskFile);

  /**
  \brief Combines near and far drawn points in one long vector to be used for edema segmentation based on region growing
  */
  VectorVectorDouble FormulateDrawingPointsForEdemaSegmentation();

  /**
  \brief Puts initial seed points in one vector to be used for tumor segmentation
  */
  VectorVectorDouble FormulateDrawingPointsForTumorSegmentation();

  /**
  \brief Save the current selected Nifti image
  */
  void SaveImage();

  //! Save an image by passing the index of the input image that needs to be saved and the appropriate file name as a QString
  void SaveImage_withFile(int indexOfInputImageToWrite, QString saveFileName);

  /**
  \brief Save the current selected DICOM image
  */
  void SaveDicomImage();

  /**
  \brief Get indices of particular label
  */
  VectorVectorDouble GetMaskLabelIndices(const int label);

  /**
  \brief Read ROI from a file into the memory
  */
  void readMaskFile(const std::string &maskFileName);

  /**
  \brief Change the size of drawing/erasing brush
  */
  void ChangeBrushSize(int size);

  /**
  \brief Changes the label in the drawing
  */
  void ChangeDrawingLabel(int drawingLabel); // multiLabel uncomment this function

  /**
  \brief Changes the selected opacity to draw with, pulling directly from the opacity selected in the drawing panel
  */
  void ChangeMaskOpacity(); // multiLabel uncomment this function

  /**
  \brief Changes the selected opacity to display by input of a number (float).
  */
  void ChangeMaskOpacity(float newOpacity);

  /**
  \brief Checks whether required images are present for the recurrence estimation application
  */
  bool CheckCompletenessOfInputData(bool & convDataPresent, bool & perfusionDataPresent, bool & dtiDataPresent, bool existingmodel);

  /**
  \brief Checks whether required images are present for the EGFRvIII estimation application
  */
  bool CheckCompletenessOfInputDataForEGFR(bool &t1ceDataPresent, bool &t2flairDataPresent, bool &perfusionDataPresent);

  /**
  \brief Sets current tissue type to be drawn
  */
  void SetActiveLandmarksType(int type, int row, int col);

  /**
  \brief on change of panel
  */
  void panelChanged(int current);

  //! Propagate the Slicer position in accordance with the selected Slicer and image IDs
  void propogateSlicerPosition(int slicerId = 0, int imageId = -1);


  /**
  \brief Sets the current selected image
  \param id The id of the current selected image
  */
  void CurrentImageChanged(std::string &id);

  /**
  \brief Updates the information in the bottom panel based on the current selected image
  */
  void ImageInfoChanged();

  /**
  \brief Close the current selected image
  */
  void CloseImage();

  /**
  \brief Close all loaded images
  */
  void CloseAllImages();

  /**
  \brief Reset the number of points in the table when all the images are closed
  */
  void ResetNumberOfPoints();

  /**
  \brief This function deals with undo
  */
  void UndoFunctionality();

  /*
  \brief Fill button logic
  */
  void FillLabel(int label);

  /**
  \brief Sets the opacity of the mask for the displayed images
  */
  void SetOpacity();

  /**
  \brief This function deals with the functionality associated with the change in the overlay slider
  */
  void overlaySliderChanged(int value);

  /**
  \brief This function deals with the functionality associated with the change in the overlay slider
  */
  void imageModalityChanged(int value);

  /**
  \brief This function deals with the functionality associated with the change in the image slider for the perfusion image
  */
  void imageSliderChanged();

  /**
  \brief Resets the transfrmation of each slicer to identity, resets the camera and renders the slicer again
  */
  void ResetTransformationToIdentity();

  /**
  \brief Renders all the images again
  */
  void UpdateRenderWindows();

  /**
  \brief Enable or disable overlay
  \param State of the overlay to set (true/false)
  */
  void overlayUseStateChanged(int state);

  /**
  \brief Passes the item, whose overlay needs to be update, to overlayChanged(QTableWidgetItem *item) function
  */
  void overlayChanged();

  /**
  \brief Sets the item selected in overlay table as overlay on item selected in images table
  */
  void overlayChanged(QTableWidgetItem *item);

  /**
  \brief Maintains record of currentlty displayed image
  */
  void CurrentPickedImageChanged(std::string id);

  /**
  \brief Displays the image currently selected in images table. Internally calls DisplayChanged(QTableWidgetItem *item) function
  */
  void DisplayChanged();

  /**
  \brief Displays the image (item) currently selected in images table
  */
  void DisplayChanged(QTableWidgetItem *item);

  /**
  \brief Updates co-ordinates based on each change in the mouse position
  \param visibity This param shows whether coordinates need to be displayed or not
  \param x LPS/RAS co-ordinate X
  \param y LPS/RAS co-ordinate Y
  \param z LPS/RAS co-ordinate Z
  \param X pixel co-ordiante X
  \param Y pixel co-ordiante Y
  \param Z pixel co-ordiante Z
  \param value Intensity value at current voxel
  */
  void MousePositionChanged(int visibility, double x, double y, double z, double X, double Y, double Z, double value);

  /**
  \brief Sets the window and level values. Calls UpdateWindowLevel function internally
  */
  void WindowLevelChanged();

  /**
  \brief Sets the window and level values. Calls UpdateWindowLevel function internally
  */
  void WindowLevelEdited();
  /**
  \brief Sets the window and level values. Calls UpdateWindowLevel function internally
  */
  void SetWindowLevel(double w, double l);
  /**
  \brief Applies the window and level values selected from spin boxes on currently displayed images
  */
  void UpdateWindowLevel();

  /**
  \brief Applies the value of threshold slider on the displayed images
  */
  void thresholdSpinBoxChanged();

  /**
  \brief Enable mask thresholding using the radio button
  */
  void EnableThresholdOfMask();

  /**
  \brief Creates a link between two currently displayed images
  \param image1 ID of the first image
  \param image2 ID of the second image
  */
  void AddLink(const std::string &image1, const std::string &image2);

  /**
  \brief Removes link between two images
  \param image1 ID of the first image
  \param image2 ID of the second image
  */
  void RemoveLink(const std::string &image1, const std::string &image2);
  /**
  \brief Sets the value of the axial, coronal or sagital sliders
  \param slicer Number of the slicer
  \param slice Value to bet set on the slider
  */
  void UpdateSlice(int slicer, int slice);
  /**
  \brief Sets the minimum and maximum values of axial, coronal or sagital sliders
  \param slicer Number of the slicer
  \param min Minimum value of the slider
  \param max Maximum value of the slider
  */
  void UpdateSliceRange(int slicer, int min, int max);

  /**
  \brief Moves the cursor on the given co-ordinates
  */
  void MoveSlicerCursor(double x, double y, double z, int mode = 0);

  std::vector<std::map<CAPTK::ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectoryForSurvival(const std::string directoryname);
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectoryForPCA(const std::string directoryname);
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectoryForPseudoProgression(const CAPTK::MachineLearningApplicationSubtype type, const std::string &directoryname, const bool &useConventionalData, const bool &useDTIData, const bool &usePerfData, const bool &useDistData);
  std::vector<std::map<CAPTK::ImageModalityType, std::string>> LoadQualifiedSubjectsFromGivenDirectoryForRecurrence(const CAPTK::MachineLearningApplicationSubtype type, const std::string &directoryname, const bool &useT1Data, const bool &useDTIData, const bool &usePerfData, const bool &useDistData);
  /**
  \brief Displays the next image in order based on keyboard input (1,2,3,...)
  */
  void ChangeImageWithOrder(SlicerManager *sm, int order);

  /**
  \brief Called internally within DisplayChange function to update link between images
  */
  void UpdateLinkManager(std::string id, int slicer, double x, double y, double z);

  /**
  \brief Called internally within DisplayChange function to update link between images
  */
  void UpdateLinkedNavigation(Slicer* refSlicer);

  //! Returns the active tab from the tab widget
  int getActiveTabId()
  {
    return m_tabWidget->currentIndex();
  }

  //! shows if the tumor points are being initialized
  bool tumorPointSelected()
  {
    return tumorPanel->mTumorPointsSelected;
  }

  //! GUI control for breast segmentation
  void ApplicationBreastSegmentation();

  //! GUI control for breast texture feature extraction
  void ApplicationTexturePipeline();

  //! GUI control for LIBRA single image mode
  void ApplicationLIBRASingle();

  //! GUI control for LIBRA batch image mode
  void ApplicationLIBRABatch();

  //! GUI control for Confetti
  void ApplicationConfetti();

  //! GUI control for SBRT Lung Field segmentation
  void ApplicationSBRTLungField();

  //! GUI control for SBRT Lung Nodule segmentation
  void ApplicationSBRTNodule();

  //! GUI control for SBRT Lung Analysis
  void ApplicationSBRTAnalysis();

  //! Convert 2D image to 3D image with a single slice and write to temp folder
  std::string ConversionFrom2Dto3D(const std::string &fileName);

  //! GUI control for Directionality Analysis
  void ApplicationDirectionality();
#ifdef BUILD_FETALBRAIN
  //! GUI control for Fetal Brain
  void ApplicationFetalBrain();
#endif

#ifdef BUILD_EGFRvIII
  //! GUI control for EGFRvIII PHI
  void ApplicationEGFR();
#endif
#ifdef BUILD_RECURRENCE
  //! GUI control for Recurrence
  void ApplicationRecurrence();
#endif
#ifdef BUILD_PSEUDOPROGRESSION
  //! GUI control for Pseudo Progression
  void ApplicationPseudoProgression();
#endif
#ifdef BUILD_ATLAS
  //! GUI control for Population Atlas
  void ApplicationPopulationAtlas();
#endif
#ifdef BUILD_ISUBTYPE
  //! GUI control for Imaging Subtype
  void ApplicationImagingSubtype();
#endif
#ifdef BUILD_MSUBTYPE
  //! GUI control for Molecular Subtype
  void ApplicationMolecularSubtype();
#endif
#ifdef BUILD_SURVIVAL
  //! GUI control for Survival
  void ApplicationSurvival();
#endif
#ifdef BUILD_EGFRvIIISVM
  //! GUI control for EGFRvIII SVM
  void ApplicationEGFRvIIISVM();
#endif
#ifdef BUILD_GEODESIC
  //! GUI control for Geodesic Segmentation
  void ApplicationGeodesic();
#endif
#ifdef BUILD_GEODESICTRAINING
  //! GUI control for Geodesic Training
  void ApplicationGeodesicTraining();
#endif
  //! GUI control for Geodesic Threshold - TBD
  void ApplicationGeodesicTreshold();
#ifdef BUILD_ITKSNAP
  //! GUI control for ITKSNAP
  void ApplicationITKSNAP();
#endif
#ifdef BUILD_WHITESTRIPE
  //! GUI control for WhiteStripe
  void ApplicationWhiteStripe();
#endif

  //! GUI control for all DeepMedic segmentations
  void ApplicationDeepMedicSegmentation(int type);

  //! GUI control for Theia
  void ApplicationTheia();

  //! GUI control for all PCA analysis
  void ApplicationPCA();
  
  //! GUI control for Perfusion Measures calculation
  void ApplicationPerfusionMeasuresCalculation();
  
  //! GUI control for Perfusion Alignment calculation
  void ApplicationPerfusionAlignmentCalculation();
  
  //! GUI control for Diffusion Measures calculation
  void ApplicationDiffusionMeasuresCalculation();

  //! GUI control for Training Module
  void ApplicationTrainingModule();

  //! Preprocessing for denoising
  void ImageDenoising();

  //! Preprocessing for mammogram preprocessing
  void ImageMamogramPreprocess();

  //! BraTS Pipeline
  void ImageBraTSPipeline();

  //! Preprocessing for bias correction
  void ImageBiasCorrection();

  //! Preprocessing for image registration
  void ImageRegistration();

  //! Preprocessing for histogram matching
  void ImageHistogramMatching();

  //! Preprocessing for Z-Scoring normalizer
  void ImageDeepMedicNormalizer();

  //! Preprocessing for skull-stripping
  void ImageSkullStripping();

  //! Preprocessing for DCM-to-NIfTI conversion
  void DCM2NIfTIConversion();

  //! Preprocessing for customized preprocessing
  void CustomPreprocessing();

  //! Enable/Disable comparison mode
  void EnableComparisonMode(bool);

  //! Get Comparison Viewers
  std::vector<vtkSmartPointer<Slicer>> GetComparisonViewers();

  //! Geodesic Training Finished Handler
  void GeodesicTrainingFinishedHandler();

  //! Geodesic Training Finished with Error Handler
  void GeodesicTrainingFinishedWithErrorHandler(QString errorMessage);

  //! Performs the registration
  void Registration(std::string fixedFileName, std::vector<std::string> inputFileNames,
    std::vector<std::string> outputFileNames, std::vector<std::string> matrixFileNames, 
    std::string metrics, bool rigidMode, bool affineMode, bool deformMode, 
    std::string radii, std::string iterations, std::string degreesOfFreedom);

  //confirm before exit
  void closeEvent(QCloseEvent * event);

  // Progress Update
  void updateProgress(int progress, std::string message = "", int max = 100);

  //! Enables "advanced mode" - no image checks are done - disabled by default
  void EnableAdvancedVisualizer()
  {
    m_advancedVisualizer = true;
  }

public:

  std::string currentPlatform;

  std::string m_tempFolderLocation;
  std::string dicomfilename; // contains the first image in the DICOM series
  std::vector<SlicerManager*> mSlicerManagers;
  std::vector<SimpleImageManager*> mNonViewingImageManager;

  QString mInputPathName;
  std::vector<QSlider*> verticalSliders;
  //
  std::string mCurrentSelectedImageId;
  std::string mCurrentPickedImageId;
  std::string mProjectVariant;
  int mCurrentPickedImageIndex;

  bool m_advancedVisualizer = false;

  //progress Bar
  QProgressBar *m_progressBar;
  QLabel* m_messageLabel;

  //
  Landmarks* mLandmarks;
  Landmarks* mSeedPoints;
  Landmarks* mTissuePoints;
  int mCurrentLandmarkTissueType;
  vtkSmartPointer<vtkImageData> mMask;

  SHAPE_MODE m_drawShapeMode;
  int mCurrentNearPoints;
  int mCurrentFarPoints;
  int mCurrentInitPoints;


  int mBorderStartX;
  int mBorderStartY;
  int mBorderStartZ;

  int mBorderEndX;
  int mBorderEndY;
  int mBorderEndZ;

  PreprocessingPipelineClass mPreprocessingObj;
  OutputWritingManager mOutputManager;
  NiftiDataManager mNifiDataManager;

  // Applications that need to be declared on fMainWindow level (multiple UI level interactions needed)
  RecurrenceEstimator mRecurrenceEstimator;
  PseudoProgressionEstimator mPseudoEstimator;
  SurvivalPredictor mSurvivalPredictor;
  EGFRvIIIIndexPredictor mEGFRvIIIPredictor;
  PopulationAtlases mPopulationAtlas;
  ImagingSubtypePredictor mImagingSubtype;
  MolecularSubtypePredictor mMolecularSubtype;
  Fetalbrain mfetalbrain;

  fHelpDialog* mHelpDlg;
  fHelpTutorial mHelpTutorial;

  std::string t1cePath;
  std::string m_imagetype_string;
  Ui::fMainWindow	m_mainWindowUI;

  std::vector<int> mActionIds;
  std::vector<int> mActionSequenceIds;

  std::vector< std::vector<PointVal> > mActionPoints;

  GenericImage::Pointer mCustomImageToThreshold;

  int mSequenceNumber, mCustomImageToThreshold_min, mCustomImageToThreshold_max;

private:
  ImageTypeFloat3D::Pointer m_InputGeomasks;
  ImageTypeShort3D::Pointer m_imgGeodesicOut;
  ImageTypeShort3D::Pointer m_imgGeodesicOutPositive;
  ImageTypeShort3D::Pointer m_imgGeodesicOutNegative;
  std::map<std::string, float> m_fetalbrainfeatures;
  int m_fetalslice;
  bool m_ComparisonMode; //! comparison mode
  vtkSmartPointer<Slicer> m_ComparisonViewerLeft, m_ComparisonViewerCenter, m_ComparisonViewerRight;

  // GeodesicTraining private variables
  GeodesicTrainingCaPTkApp<2>* m_GeodesicTrainingCaPTkApp2D;
  GeodesicTrainingCaPTkApp<3>* m_GeodesicTrainingCaPTkApp3D;
  std::string m_GeodesicTrainingFirstFileNameFromLastExec = "";
  bool m_IsGeodesicTrainingRunning = false;

  std::thread m_ExternalProcessThread;

  struct DicomDictTagAndVal
  {
    std::string tag;
    std::string value;

    DicomDictTagAndVal(const std::string &input_tag, const std::string &input_value) :
      tag(input_tag), value(input_value)
    { }

    DicomDictTagAndVal(const std::string &input_tag) :
      tag(input_tag)
    { }

    void SetValue(const std::string &input_value)
    {
      value = input_value;
    }
  };

  void RegistrationWorker(std::vector<std::string> compVector, std::vector<std::string> inputFileNames,
    std::vector<std::string> outputFileNames, std::vector<std::string> matrixFileNames);

  bool m_skipTutorialOnNextRun = false;

  std::atomic<int> m_NumberOfUnfinishedExternalProcesses = { 0 };

  int  startExternalProcess(const QString &application, const QStringList &arguments);
};
#if __GNUC__
#pragma GCC visibility pop
#endif
//-------------------------------------------------------------------------------------

#endif
