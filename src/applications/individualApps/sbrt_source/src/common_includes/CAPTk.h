/**
\file CAPTk.h

\brief Main header location for CAPTk
*/

#pragma once

// system includes
#include <algorithm>
#include <errno.h>
#include <locale.h>
#include <ios>
#include <fstream>
#include <sstream>
#include <iostream>
#include <ostream>
#include <utility>
#include <cassert>
#include <ctime>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

// ITK
#include "itkImage.h"
#include "itkVariableLengthVector.h"
#include "itkVariableSizeMatrix.h"
#include "itkFlipImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkDiffusionTensor3DReconstructionImageFilter.h"
#include "itkTensorFractionalAnisotropyImageFilter.h"
#include "itkTensorRelativeAnisotropyImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
//#include "itkImageToVectorImageFilter.h" // This class is deprecated. You should use itkComposeImageFilter instead.
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageMaskSpatialObject.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkByteSwapper.h"
#include "itkCommand.h"
#include "itkNumericSeriesFileNames.h"
#include "itkCastImageFilter.h"
#include "itkVTKImageImport.h"
#include "itkImageDuplicator.h"
//#include "itkNiftiImageIO.h"
#include "itkImageRegionIterator.h"
#include "itkMetaDataObject.h"
#include "itkImageSeriesWriter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#if ITK_VERSION_MAJOR >= 4
#include "itkMultiplyImageFilter.h"
#else
#include "itkMultiplyByConstantImageFilter.h"
#endif
#if (ITK_VERSION_MAJOR < 4) || defined(ITKV3_COMPATIBILITY)
#include "itkAnalyzeImageIO.h"
#endif

// VTK
#include "vtksys/SystemTools.hxx"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkTransform.h"
#include "vtkInteractorStyle.h"
#include "vtkRenderWindow.h"
#include "vtkPolyData.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPropPicker.h"
#include "vtkRenderer.h"
#include "vtkImageActor.h"
#include "vtkCommand.h"
#include "vtkTextProperty.h"
#include "vtkTextActor.h"
#include "vtkTextSource.h"
#include "vtkActor2D.h"
#include "vtkCursor2D.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkProperty2D.h"
#include "vtkCornerAnnotation.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkToolkits.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkClipPolyData.h"
#include "vtkGlyph3D.h"
#include "vtkMath.h"
#include "vtkCursor3D.h"
#include "vtkProperty.h"
#include "vtkLight.h"
#include "vtkLightCollection.h"
#include "vtkScalarBarActor.h"
#include "vtkLookupTable.h"
#include "vtkLabeledDataMapper.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkRendererCollection.h"
#include "vtkCamera.h"
#include "vtkCallbackCommand.h"
#include "vtkPolyDataMapper.h"
#include "vtkBox.h"
#include "vtkExtractVOI.h"
#include "vtkSphereSource.h"
#include "vtkCutter.h"
#include "vtkAssignAttribute.h"
#include "vtkImageAccumulate.h"
#include "vtkImageReslice.h"
#include "vtkCubeSource.h"
#include "vtkImageViewer2.h"
#include "vtkImageMapToColors.h"
#include "vtkWindowToImageFilter.h"
#include "vtkMatrix4x4.h"
#include "vtkRegularPolygonSource.h"
#include "vtkImageExport.h"
#include "vtkDoubleArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPCAStatistics.h"
#include "vtkStringArray.h"
#include "vtkTable.h"
#if VTK_MAJOR_VERSION >= 6 || (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 10)
#include "vtkImageMapper3D.h"
#endif
#if VTK_MAJOR_VERSION >= 6
#include "vtkInformation.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkAlgorithmOutput.h"
#endif

// Qt4
#include <QtGlobal>
#include <QApplication>
#include <QPixmap>
#include <QSplashScreen>
#include <QTimer>
#include <QDesktopWidget>
#include <QDir>
#include <QObject>
#include <QFileDialog>
#include <QMessageBox>
#include <QInputDialog>
#include <QUrl>
#include <QSettings>
#include <QEvent>
#include <QTreeWidget>
#include <QTableWidget>
#include <QProgressDialog>
#include <QPushButton>
#include <QLineEdit>
#include <QCheckBox>
#include <QSpinBox>
#include <QProcess>
#include <QListView>
#include <QDebug>
#include <QDirIterator>
#include <QDragEnterEvent>
#include <QDir>
#include "Qt/qcombobox.h"
#include "Qt/qrect.h"
#include "qmessagebox.h"

#ifdef Q_OS_WIN
#include <windows.h> // for Sleep
#endif

#include "cbicaUtilities.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"

//const int DT_UNKNOWN = 0;


#define ROUND(t)			((t) < 0 ? (int) ((t)-0.5f) : (int) ((t)+0.5f))
#define CLAMP(a, l, h)		((a) < (l) ? (l) : (a) > (h) ? (h) : (a))
#ifndef __max
#define __max(a,b)			(((a) > (b)) ? (a) : (b))
#endif
#ifndef __min
#define __min(a,b)			(((a) < (b)) ? (a) : (b))
#endif


// can all the following be replaced with enums (which would make the code significantly clearer)? Some examples are given below...

// let us use a single format, and ensure proper conversion (with messages telling users what is going on) between one and the other
//#define USE_LANDMARKS_LPS_COORD
#define USE_LANDMARKS_RAS_COORD // [TBD] - this isn't needed since we will be working with a single coordinate system in the future

#define SVM_MODEL_EXTENSIONS "Models (*.model)"

#define PERF_POINTS_ESTIM 10 // [TBD] - should this be hard-coded

//enum ImageExtension
//{
//  NIFTI, DICOM
//};
#define IMAGE_NIFTI 0 // [TBD] - convert to enum
#define IMAGE_DICOM 1

//enum ImageModalityType
//{
//  UNDEFINED, T1, T1CE, T2, T2FLAIR, PERFUSION, DTI, RECURRENCE_OUTPUT
//};
#define IMAGE_TYPE_UNDEFINED			0 // [TBD] - convert to enum
#define IMAGE_TYPE_T1					1
#define IMAGE_TYPE_T1CE					2
#define IMAGE_TYPE_T2					3
#define IMAGE_TYPE_T2FLAIR				4
#define IMAGE_TYPE_PERFUSION			5
#define IMAGE_TYPE_DTI					6
#define IMAGE_TYPE_RECURRENCE_OUTPUT	7
#define IMAGE_TYPE_PP			8


/enum NV_Image
//{
//  T1_ROW, T2_ROW, PERF_ROW
//};
#define NV_IMAGE_T1_ROW 0 // [TBD] - convert to enum
#define NV_IMAGE_T2_ROW 1
#define NV_IMAGE_PERFUSION_ROW 2

//enum MaskProcess
//{
//  RECURRENCE_ROW, NON_RECURRENCE_ROW, TOTAL_ROW
//};
#define MASK_RECURRENCE_ROW 0 // [TBD] - convert to enum
#define MASK_NON_RECURRENCE_ROW 1
#define MASK_TOTAL_ROW 2

//enum PreProcessType
//{
//  SUSAN, BIAS_CORRECT, REGISTRATION, RECURRENCE_MAP
//};
#define WRITING_TYPE_SUSAN 0 // [TBD] - convert to enum
#define WRITING_TYPE_BIASCORRECTION 1
#define WRITING_TYPE_REGISTRATION 2
#define WRITING_TYPE_RECURRENCEMAP 3

//enum TrainingLabel
//{
//  NEAR, FAR
//};
#define TRAINING_LABEL_NEAR 0 // [TBD] - convert to enum
#define TRAINING_LABEL_FAR 1

//enum RunningMode
//{
//  NORMAL, CLUSTER
//};
#define RUNNING_MODE_NORMAL 1 // [TBD] - convert to enum
#define RUNNING_MODE_CLUSTER 2


//enum LabelPointType
//{
//  NEAR_POINT, FAR_POINT, INIT_POINT
//};
#define NEAR_POINT_LABEL 1 // [TBD] - convert to enum
#define FAR_POINT_LABEL 2
#define INIT_POINT_LABEL 3

//enum PixelValuesForROI
//{
//  INIT_POINT = 100, NEAR_POINT = 150, FAR_POINT = 255
//};
#define NEAR_POINT_SAVE_LABEL 150 // why is this and INIT the same value?
#define INIT_POINT_SAVE_LABEL 150 // [TBD] - convert to enum
#define FAR_POINT_SAVE_LABEL 255

#define LANDMARK_MAX_NearFarDrawingManager 200 // what does this do? // [TBD] - convert to enum
#define MAX_SEED_DRAWING_POINTS 200 // what does this do?
#define MAX_NEARFAR_REGION_POINTS 1000 // what does this do?

#define NO_OF_PCS 5 // what does this do? // [TBD] - convert to enum
#define NO_OF_FEATURES 13 // what does this do?
#define NO_OF_PCA_FEATURES 45 // what does this do?

#define EGFR_PCS 3 // what does this do?

// Valid Image extensions
#define IMAGES_EXTENSIONS "Images (*.nii.gz *.nii *dcm)" // why this and not the ones below?
#define NII_EXT ".nii" // [TBD] - convert to enum
#define NII_GZ_EXT ".nii.gz"
#define HDR_EXT ".hdr"
#define IMG_EXT ".img"


// Common defines - these are required for GCC compilation 
typedef itk::VariableLengthVector< double > VariableLengthVectorType;
typedef itk::VariableSizeMatrix< double > VariableSizeMatrixType;
typedef std::vector< std::vector < double > > VectorVectorDouble;
typedef std::vector < double > VectorDouble;
typedef itk::Image< float, 3 > ImageTypeFloat3D;
typedef itk::Image< float, 4 > ImageTypeFloat4D;
typedef itk::Image< short, 3 > ImageTypeShort3D;


enum DrawingActions
{
  NEAR_DRAWING = 1,
  NEAR_ERASING = 2,
  FAR_DRAWING = 3,
  FAR_ERASING = 4
};

struct NonNativeApp
{
  std::string name;
  std::string path;

  //! Default Constructor
  NonNativeApp()
  {
    name = "";
    path = "";
  }

  //! Default Constructor with values
  NonNativeApp(const std::string &inputName, const std::string &inputPath) :
    name(inputName), path(inputPath){};
};

const std::string loggerFolder = cbica::getUserHomeDirectory() + "/CaPTk-" + std::string(PROJECT_VERSION) + "_logs/";
const std::string loggerFile = loggerFolder + cbica::replaceString(cbica::getCurrentLocalDate(), ":", "-") + ".log";
const std::string aboutScreenSeen = loggerFolder + "aboutSeen.txt"; // file to check if the user has seen "about" page or not

static QRect screenSize;

//! Initialize the screen width of the current machine
inline void SetScreenSize(const QRect &inputScreenSize)
{
  screenSize = inputScreenSize;
}

//! Cross platform Sleep
inline void qSleep(int ms)
{
  if (ms <= 0)
  {
    cbica::Logging(loggerFile, "Sleep time of zero or less is not defined.");
    exit(EXIT_FAILURE);
  }

#ifdef Q_OS_WIN
  Sleep(uint(ms));
#else
  struct timespec ts = { ms / 1000, (ms % 1000) * 1000 * 1000 };
  nanosleep(&ts, NULL);
#endif
}

template< class QThingType >
inline void fixSize(QThingType *thingToFix, QRect &screenSize)
{
  thingToFix->setMinimumWidth(100); // needs to be screenSize dependent 
  thingToFix->setMinimumHeight(50); // needs to be screenSize dependent 
  thingToFix->setMaximumWidth(std::round(static_cast<long>((screenSize.width() + 1) / 10)));
  thingToFix->setMaximumHeight(std::round(static_cast<long>((screenSize.height() + 1) / 10)));
}

inline void fixComboBox(QComboBox *comboBoxToEdit, QRect &screenSize)
{
  fixSize< QComboBox >(comboBoxToEdit, screenSize);
  //comboBoxToEdit->setEditable(true);
  //comboBoxToEdit->lineEdit()->setAlignment(Qt::AlignCenter);
  //comboBoxToEdit->lineEdit()->setReadOnly(true);
  for (int i = 0; i < comboBoxToEdit->count(); i++)
  {
    comboBoxToEdit->setItemData(i, Qt::AlignCenter, Qt::TextAlignmentRole);
  }
}

inline void fixComboBox(QComboBox *comboBoxToEdit)
{
  //comboBoxToEdit->setEditable(true);
  //comboBoxToEdit->lineEdit()->setAlignment(Qt::AlignCenter);
  //comboBoxToEdit->lineEdit()->setReadOnly(true);
  for (int i = 0; i < comboBoxToEdit->count(); i++)
  {
    comboBoxToEdit->setItemData(i, Qt::AlignCenter, Qt::TextAlignmentRole);
  }
}

/**
\brief Displays the message in a small message box with default window title as "Error"

This function is used for displaying error messages. For success, use the statusBar()->showMessage() function.

\param message The message to be displayed in the message box
\param windowTitle The title of the message box, defaults to "Error"
*/
inline void ShowErrorMessage(const std::string &message, const std::string &windowTitle = "Error")
{
  QMessageBox *box = new QMessageBox();
  box->setIcon(QMessageBox::Information);
  box->addButton(QMessageBox::Ok);
  box->setText(QString(message.c_str()));
  box->setWindowTitle(QString(windowTitle.c_str()));
  box->exec();
}