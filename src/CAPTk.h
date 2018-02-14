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
#include <math.h>
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
#include "itkChangeInformationImageFilter.h"
//#include "itkNiftiImageIO.h"
#include "itkImageRegionIterator.h"
#include "itkMetaDataObject.h"
#include "itkImageSeriesWriter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkCSVNumericObjectFileWriter.h"
#include "itkCSVArray2DFileReader.h"
#include "itkConnectedComponentAlgorithm.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
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
#include "qstring.h"

// OpenCV
#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"
#include "opencv2/imgproc.hpp"

#include <ostream>
#include <fstream>

#ifdef Q_OS_WIN
#include <windows.h> // for Sleep
#include <io.h>
#endif

#include "cbicaUtilities.h"
#include "cbicaLogging.h"
#include "cbicaITKSafeImageIO.h"

//const int DT_UNKNOWN = 0;
//#define _USE_MATH_DEFINES

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

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
//#define USE_LANDMARKS_LPS_COORD // both of these are redundant since the current saving of tissue and tumor points is happening in the real coordinate system
#define USE_LANDMARKS_RAS_COORD // [TBD] - this isn't needed since we will be working with a single coordinate system in the future

#define SVM_MODEL_EXTENSIONS "Models (*.model)"

enum ImageExtension
{
  NIfTI, DICOM
};

enum ImageModalityType
{
	IMAGE_TYPE_UNDEFINED, IMAGE_TYPE_T1, IMAGE_TYPE_T1CE, IMAGE_TYPE_T2, IMAGE_TYPE_T2FLAIR, IMAGE_TYPE_AX, IMAGE_TYPE_FA, IMAGE_TYPE_RAD, IMAGE_TYPE_TR, IMAGE_TYPE_PERFUSION, 
  IMAGE_TYPE_DTI, IMAGE_TYPE_RECURRENCE_OUTPUT, IMAGE_TYPE_PP, IMAGE_TYPE_CT, IMAGE_TYPE_PET, IMAGE_TYPE_PSR, IMAGE_TYPE_PH, 
  IMAGE_TYPE_RCBV, IMAGE_TYPE_SEG, IMAGE_TYPE_ATLAS, IMAGE_TYPE_PARAMS, IMAGE_TYPE_SUDOID, IMAGE_TYPE_NEAR, IMAGE_TYPE_FAR, IMAGE_MAMMOGRAM, IMAGE_TYPE_FEATURES
};

static const char ImageModalityString[IMAGE_TYPE_FEATURES + 1][15] = { "DEF", "T1", "T1Gd", "T2", "FLAIR", "DTI_AX", "DTI_FA", "DTI_RAD", "DTI_TR", "PERFUSION",
"DTI", "REC", "PP", "CT", "PET", "pSR", "PH", 
"RCBV", "SEG", "ATLAS", "PARAMS", "SUDOID", "NEAR", "FAR", "FFDM", "FEAT" };

enum MachineLearningApplicationSubtype
{
	TRAINING, TESTING
};

enum ApplicationCallingSVM
{
  Default, Recurrence, Survival, SBRT  
};

enum SURVIVAL_SVM_PARAMS
{
	BESTC_6MONTH = 32, BESTG_6MONTH = 1/32, BESTC_18MONTH = 1/2, BESTG_18MONTH = 1/32
};

enum GLISTR_OUTPUT_LABELS
{
	EDEMA_1 = 100, TUMOR_1 = 200, VENT_1 = 10, NONENHANCING_1 = 175, EDEMA_2 = 5,   TUMOR_2 = 7,   VENT_2 = 1,  NONENHANCING_2 = 6, ALL = 0
};

enum MOLECULAR_SUBTYPES
{
	PRONEURAL = 1, NEURAL = 2, MESSENCHYMAL = 3, CLASSICAL = 4
};

enum VOXEL_STATUS
{
	OFF, ON
};


//enum RunningMode
//{
//  NORMAL, CLUSTER
//};

// Valid Image extensions
//#define IMAGES_EXTENSIONS "Images (*.nii.gz *.nii *dcm)" // why this and not the ones below?
#define NII_EXT ".nii" // [TBD] - convert to enum
#define NII_GZ_EXT ".nii.gz"
#define HDR_EXT ".hdr"
#define IMG_EXT ".img"
#define PARAM_EXT ".txt"
#define CSV_EXT ".csv"

static QString IMAGES_EXTENSIONS = "Images (*.nii.gz *.nii *dcm)";

// Common defines - these are required for GCC compilation 
using VariableLengthVectorType = itk::VariableLengthVector< double >; // TBD: move from itk::matrix to vnl::matrix for simplicity
using VariableSizeMatrixType = itk::VariableSizeMatrix< double >; // TBD: move from itk::matrix to vnl::matrix for simplicity
using VectorVectorDouble = std::vector< std::vector < double > >;
using VectorDouble = std::vector < double >;
using ImageTypeFloat3D = itk::Image< float, 3 >;
using ImageTypeFloat4D = itk::Image< float, 4 >;
using ImageTypeShort3D = itk::Image< short, 3 >;
using ImageTypeFloat3DIterator = itk::ImageRegionIteratorWithIndex< ImageTypeFloat3D >;

//! Structure to define a point value to check if it is defined in the image or not
struct PointVal
{
	PointVal()
	{
		x = y = z = value = 0;
	}
	PointVal(int _x, int _y, int _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}
	PointVal(int _x, int _y, int _z, int _value)
	{
		x = _x;
		y = _y;
		z = _z;
		value = _value;
	}
	bool isWithinRange(const int* dims)
	{
		if (x >= 0 && x < dims[0] && y >= 0 && y < dims[1] && z >= 0 && z < dims[2])
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool isValid()
	{
		if (x >= 0 && y>= 0 && z >= 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	PointVal getInvalidPt()
	{
		return PointVal(-1, -1, -1, -1);
	}
	int x, y, z, value;
};

Q_DECLARE_METATYPE(PointVal);
Q_DECLARE_METATYPE(std::vector< PointVal >);

//! Structure to store non-native (read: standalone applications) for easier parsing
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

const std::string loggerFolderBase = cbica::getUserHomeDirectory() + "/." + std::string(PROJECT_NAME) + "/"; // kept separate because of folder creation
const std::string loggerFolder = loggerFolderBase + std::string(PROJECT_VERSION) + "/";
const std::string loggerFile = loggerFolder + cbica::replaceString(cbica::getCurrentLocalDate(), ":", "-") + ".log";
const std::string aboutScreenSeen = loggerFolder + "aboutSeen.txt"; // file to check if the user has seen "about" page or not
const std::string tutorialScreen = loggerFolder + "tutorial.txt"; // file to check if the user has seen "about" page or not
const std::string closeConfirmation = loggerFolder + "close.txt"; // file to check if the user has seen "about" page or not


inline void fixComboBox(QComboBox *comboBoxToEdit)
{
  for (int i = 0; i < comboBoxToEdit->count(); i++)
  {
    comboBoxToEdit->setItemData(i, Qt::AlignCenter, Qt::TextAlignmentRole);
  }
}

/**
\brief Displays the message in a small message box with default window title as "Error"

This function is used for displaying error messages. 

For detailed output results, use ShowMessage() function. For simple success notification, use the updateProgress function.

\param message The message to be displayed in the message box
\param windowTitle The title of the message box, defaults to "Error"
*/
inline void ShowErrorMessage(const std::string &message, QWidget *boxParent = NULL, const std::string &windowTitle = "Error")
{
  if (boxParent == NULL)
  {
    QMessageBox *box = new QMessageBox();
    box->setIcon(QMessageBox::Information);
    box->addButton(QMessageBox::Ok);
    box->setText(QString(message.c_str()));
    box->setWindowTitle(QString(windowTitle.c_str()));
    box->exec();
  }
  else
  {
    QMessageBox *box = new QMessageBox(boxParent);
    box->setIcon(QMessageBox::Information);
    box->addButton(QMessageBox::Ok);
    box->setText(QString(message.c_str()));
    box->setWindowTitle(QString(windowTitle.c_str()));
    box->exec();
  }
}

/**
\brief Displays the message in a small message box with default window title as "Output Results"

This function is used for displaying regular output results which aren't errors. 

For error, use ShowErrorMessage() function. For simple success notification, use the updateProgress(), function.

\param message The message to be displayed in the message box
\param windowTitle The title of the message box, defaults to "Output Results"
*/
inline void ShowMessage(const std::string &message, const std::string &windowTitle = "Output Results")
{
  QMessageBox *box = new QMessageBox();
  box->setIcon(QMessageBox::Information);
  box->setStandardButtons(QMessageBox::Ok);
  box->setText(QString(message.c_str()));
  box->setWindowTitle(QString(windowTitle.c_str()));
  box->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
  //box->setModal(false); // if you want it non-modal
  box->setWindowModality(Qt::NonModal);
  QCoreApplication::processEvents();
  box->show();
  QCoreApplication::processEvents();
  //box->open(this, SLOT(msgBoxClosed(QAbstractButton*)));
}

/**
\brief Train and save the SVM classifier

\param trainingDataAndLabels Input training data with last column as training labels
\param outputModelName File to save the model
*/
inline VectorDouble trainOpenCVSVM(const VariableSizeMatrixType &trainingDataAndLabels, const std::string &outputModelName, bool considerWeights, int ApplicationCallingSVM)
{
  cv::Mat trainingData = cv::Mat::zeros(trainingDataAndLabels.Rows(), trainingDataAndLabels.Cols()-1, CV_32FC1), 
    trainingLabels = cv::Mat::zeros(trainingDataAndLabels.Rows(), 1, CV_32FC1);

  //// debugging purposes
  //std::ofstream file;
  //std::string base, path, ext;
  //cbica::splitFileName(outputModelName, path, base, ext);
  //file.open((path + base + "_trainingAndLabels.csv").c_str());
  //for (size_t i = 0; i < trainingDataAndLabels.Rows(); i++)
  //{
  //  for (size_t j = 0; j < trainingDataAndLabels.Cols(); j++)
  //  {
  //    file << trainingDataAndLabels(i, j) << ",";
  //  }
  //  file << "\n";
  //}
  //file.close();

  cv::Mat trainingWeights = cv::Mat::zeros(2, 1, CV_32FC1);
  size_t label1Counter = 0, label2Counter = 0;
  // fast cv::Mat access 
  for (unsigned int i = 0; i < trainingDataAndLabels.Rows(); ++i)
  {
    for (unsigned int j = 0; j < trainingDataAndLabels.Cols() - 1; ++j)
    {
      trainingData.ptr< float >(i)[j] = trainingDataAndLabels(i, j);
    }

    // last column consists of labels
    trainingLabels.ptr< float >(i)[0] = trainingDataAndLabels(i, trainingDataAndLabels.Cols() - 1);
    if (considerWeights)
    {
      // start counter to assign weights
      if (trainingLabels.ptr< float >(i)[0] == 0)
      {
        label1Counter++;
      }
      else
      {
        label2Counter++;
      }
    }
  }

  // slow cv::Mat access 
  //for (size_t i = 0; i < trainingDataAndLabels.Rows(); i++)
  //{
  //  for (size_t j = 0; j < trainingDataAndLabels.Cols() - 1; j++)
  //  {
  //    trainingData.at< float >(i, j) = trainingDataAndLabels(i, j);
  //  }
  //
  //  // last column consists of labels
  //  trainingLabels.at< float >(i,0) = trainingDataAndLabels(i, trainingDataAndLabels.Cols() - 1);
  //
  //  if (considerWeights)
  //  {
  //    // start counter to assign weights
  //    if (trainingDataAndLabels(i, trainingDataAndLabels.Cols() - 1) == 0)
  //    {
  //      label1Counter++;
  //    }
  //    else
  //    {
  //      label2Counter++;
  //    }
  //  }
  //}
  //cv::Mat diff;
  //cv::compare(trainingData, trainingData_temp, diff, cv::CMP_NE);
  //int nz = cv::countNonZero(diff);

  //cv::compare(trainingLabels, trainingLabels_temp, diff, cv::CMP_NE);
  //int nz2 = cv::countNonZero(diff);

  trainingLabels.convertTo(trainingLabels, CV_32SC1);
  
  auto svm = cv::ml::SVM::create();
  svm->setType(cv::ml::SVM::C_SVC);
  //svm->setC(1);
  //svm->setGamma(0.01);
  // parameters for auto_train
  cv::ml::ParamGrid grid_C(-5, 5, 2); // Parameter C of a SVM optimization problem (C_SVC / EPS_SVR / NU_SVR)
  cv::ml::ParamGrid grid_Gamma(-5, 5, 2); // Parameter \gamma of a kernel function (POLY / RBF / SIGMOID / CHI2)
  cv::ml::ParamGrid grid_P(-5, 5, 2); // Parameter \epsilon of a SVM optimization problem (EPS_SVR)
  cv::ml::ParamGrid grid_Nu(-5, 5, 2); // Parameter \nu of a SVM optimization problem (NU_SVC / ONE_CLASS / NU_SVR)
  cv::ml::ParamGrid grid_Degree(0, 5, 1); // Parameter degree of a kernel function (POLY)
  cv::ml::ParamGrid grid_Coeff0(-5, 5, 2); // this is Parameter coef0 of a kernel function (POLY / SIGMOID)
  if (ApplicationCallingSVM == Recurrence)
  {
    svm->setKernel(cv::ml::SVM::RBF); // using this produces terrible results for recurrence
  }
  else
  {
    svm->setKernel(cv::ml::SVM::LINEAR);
  }
  svm->setC(1);
  svm->setGamma(0.01);
  //svm->setKernel(cv::ml::SVM::POLY); // this crashes

  // assign penalties
  //if (considerWeights)
  //{
  //  trainingWeights.ptr< float >(0)[0] = label1Counter;
  //  trainingWeights.ptr< float >(1)[0] = label2Counter;
  //  //trainingWeights.at< float >(0, 0) = label1Counter;
  //  //trainingWeights.at< float >(1, 0) = label2Counter;
  //  svm->setClassWeights(trainingWeights);
  //}
  bool res = true;
  std::string msg;

  try
  {
    res = svm->train(trainingData, cv::ml::ROW_SAMPLE, trainingLabels);
    //res = svm->trainAuto(cv::ml::TrainData::create(trainingData, cv::ml::ROW_SAMPLE, trainingLabels), 
     // 10, grid_C, grid_Gamma, grid_P, grid_Nu, grid_Coeff0, grid_Degree, true);
  }
  catch (cv::Exception ex)
  {
     msg = ex.what();
  }
  //---------------------------------find distances on training adat----------------------------------
  VectorDouble returnVec;
  returnVec.resize(trainingData.rows);
  cv::Mat predicted(1, 1, CV_32F); // not sure if this is required if we do train_auto for cross-validation
  if (ApplicationCallingSVM == Recurrence) // only recurrence needs the distance map
  {
    for (int i = 0; i < trainingData.rows; i++)
    {
      cv::Mat sample = trainingData.row(i);
      //float value = svm->predict(sample, cv::Mat(), cv::ml::StatModel::RAW_OUTPUT);`
      /*float value = */svm->predict(sample, predicted, true /*cv::ml::StatModel::RAW_OUTPUT*/);
      returnVec[i] = /*predicted.at<float>(0, 0)*/predicted.ptr< float >(0)[0];
    }
  }
  else
  {
    returnVec.push_back(1.0); // just so survival application doesn't thrown an error
  }
  //--------------------------------------------------------------------------------------------------
  //svm->setTermCriteria(cv::TermCriteria(cv::TermCriteria::MAX_ITER, (int)1e7, 1e-6));
  if (res)
  {
    //cv::Mat sv = svm->getUncompressedSupportVectors();
    svm->save(outputModelName);
  }
  else
  {
    //return false;
  }

  return returnVec;
}

/**
\brief Load the SVM classifier and get the distances from the hyperplane

\param testingData Input training data with last column as training labels
\param inputModelName File to save the model
\return Distances of classification
*/
inline VectorDouble testOpenCVSVM(const VariableSizeMatrixType &testingData, const std::string &inputModelName)
{
  auto svm = cv::Algorithm::load<cv::ml::SVM>(inputModelName);

  //std::ofstream file;
  //file.open("Z:/Projects/testingData.csv");
  //for (size_t i = 0; i < testingData.Rows(); i++)
  //{
  //  for (size_t j = 0; j < testingData.Cols(); j++)
  //  {
  //    file << testingData(i, j) << ",";
  //  }
  //  file << "\n";
  //}
  //file.close();

  VectorDouble returnVec;
  //VariableSizeMatrixType returnMat;
  //returnMat.SetSize(testingData.Rows(), 1);
  //returnVec.resize(testingData.Rows());
  cv::Mat testingDataMat = cv::Mat::zeros(testingData.Rows(), testingData.Cols() - 1, CV_32FC1), outputProbs;
  
  // fast cv::Mat access
  for (unsigned int i = 0; i < testingData.Rows(); ++i)
  {
    for (unsigned int j = 0; j < testingData.Cols(); ++j)
    {
      testingDataMat.ptr< float >(i)[j] = testingData(i, j);
    }
  }

  // slow cv::Mat access
  //for (size_t i = 0; i < testingData.Rows(); i++)
  //{
  //  for (size_t j = 0; j < testingData.Cols() - 1; j++)
  //  {
  //    testingDataMat.at< float >(i, j) = testingData(i, j);
  //  }
  //}
  
  // see http://docs.opencv.org/trunk/db/d7d/classcv_1_1ml_1_1StatModel.html#af1ea864e1c19796e6264ebb3950c0b9a for details regarding why '1'
  //svm->predict(testingDataMat, returnVec, 1);

  returnVec.resize(testingDataMat.rows);
  cv::Mat predicted(1, 1, CV_32F);
  for (int i = 0; i < testingDataMat.rows; i++)
  {
    cv::Mat sample = testingDataMat.row(i);
    /*float value = */svm->predict(sample, predicted, true/*cv::ml::StatModel::RAW_OUTPUT*/);
    //float p = /*predicted.at<float>(0, 0)*/predicted.ptr< float >(0)[0];
    returnVec[i] = predicted.ptr< float >(0)[0];
  }

  //for (size_t i = 0; i < outputProbs.rows; i++)
  //{
  //  returnMat[i] = outputProbs.at<float>(i, 0);
  //}

  return returnVec;
}

/**
\brief Load the SVM classifier and predict the probabilities

\param testingData Input training data with last column as training labels
\param inputModelName File to save the model
\return Probabilities of classification
*/
//inline VectorDouble testOpenCVSVM_Probs(const VariableSizeMatrixType &testingData, const std::string &inputModelName)
//{
//  auto svm = cv::Algorithm::load<cv::ml::SVM>(inputModelName);
//
//  VectorDouble returnVec;
//  returnVec.resize(testingData.Rows());
//  cv::Mat testingDataMat = cv::Mat::zeros(testingData.Rows(), testingData.Cols(), CV_32FC1), outputProbs;
//
//  //// fast cv::Mat access -- needs to be checked
//  //float *p;
//  //for (int i = 0; i < testingData.Rows(); ++i)
//  //{
//  //  p = testingDataMat.ptr< float >(i);
//  //  for (int j = 0; j < testingData.Cols(); ++j)
//  //  {
//  //    p[j] = testingData(i, j);
//  //  }
//  //}
//
//  for (size_t i = 0; i < testingData.Rows(); i++)
//  {
//    for (size_t j = 0; j < testingData.Cols(); j++)
//    {
//      testingDataMat.at< float >(i, j) = testingData(i, j);
//    }
//  }
//
//  // see http://docs.opencv.org/trunk/db/d7d/classcv_1_1ml_1_1StatModel.html#af1ea864e1c19796e6264ebb3950c0b9a for details regarding why '1'
//  svm->predict(testingDataMat, returnVec, 1);
//
//  float max_1, max_2;
//  for (size_t i = 0; i < testingDataMat.rows; i++)
//  {
//    float temp;
//    //svm->
//  }
//
//  for (size_t i = 0; i < returnVec.size(); i++)
//  {
//    //returnVec[i] = 1.0 / (1.0 + exp( returnVec[i])); // gives probs for class 1
//    returnVec[i] = 1.0 / (1.0 + exp(-returnVec[i])); // gives probs for class 2
//  }
//  
//  return returnVec;
//}

/**
\brief Guess Image Type

\param str String to guess
\return deduced type
*/
inline int guessImageType(const std::string &fileName)
{
  int ImageSubType = IMAGE_TYPE_UNDEFINED;
  std::string fileName_wrap = fileName;
  std::transform(fileName_wrap.begin(), fileName_wrap.end(), fileName_wrap.begin(), ::tolower);
  if ((fileName_wrap.find("_t1ce") != std::string::npos) || (fileName_wrap.find("_t1-gad") != std::string::npos) ||
    (fileName_wrap.find("_t1-ce") != std::string::npos) || (fileName_wrap.find("_t1-gd") != std::string::npos) ||
    (fileName_wrap.find("_t1gd") != std::string::npos) || (fileName_wrap.find("t1gd_") != std::string::npos) ||
    (fileName_wrap.find("t1ce") != std::string::npos) || (fileName_wrap.find("t1-gad") != std::string::npos) ||
    (fileName_wrap.find("t1-ce") != std::string::npos) || (fileName_wrap.find("t1-gd") != std::string::npos))
  {
    ImageSubType = IMAGE_TYPE_T1CE;
  }
  else if ((fileName_wrap.find("_t1") != std::string::npos) || (fileName_wrap.find("t1_") != std::string::npos))
  {
    ImageSubType = IMAGE_TYPE_T1;
  }
  else if ((fileName_wrap.find("_t2") != std::string::npos) || (fileName_wrap.find("t2_") != std::string::npos))
  {
    if ((fileName_wrap.find("flair") != std::string::npos))
    {
      ImageSubType = IMAGE_TYPE_T2FLAIR;
    }
    else
    {
      ImageSubType = IMAGE_TYPE_T2;
    }
  }
  else if ((fileName_wrap.find("_flair") != std::string::npos) || (fileName_wrap.find("flair_") != std::string::npos))
  {
    ImageSubType = IMAGE_TYPE_T2FLAIR;
  }
  else if ((fileName_wrap.find("_dti") != std::string::npos) || (fileName_wrap.find("dti_") != std::string::npos) ||
    (fileName_wrap.find("_ad") != std::string::npos) || (fileName_wrap.find("ad_") != std::string::npos) || // dti sub-modalities
    (fileName_wrap.find("_ax") != std::string::npos) || (fileName_wrap.find("ax_") != std::string::npos) ||
    (fileName_wrap.find("_b0") != std::string::npos) || (fileName_wrap.find("b0_") != std::string::npos) ||
    (fileName_wrap.find("_fa") != std::string::npos) || (fileName_wrap.find("fa_") != std::string::npos) ||
    (fileName_wrap.find("_tr") != std::string::npos) || (fileName_wrap.find("tr_") != std::string::npos))
  {
    ImageSubType = IMAGE_TYPE_DTI;
  }
  else if ((fileName_wrap.find("_ct2pet.") != std::string::npos) || (fileName_wrap.find("_ct.") != std::string::npos))
  {
    ImageSubType = IMAGE_TYPE_CT;
  }
  else if ((fileName_wrap.find("_pet.") != std::string::npos) )
  {
    ImageSubType = IMAGE_TYPE_PET;
  }
  return ImageSubType;
}

/**
\brief Open existing directory (this is the mechanism to use for saving directories as well)

\param parent The QWidget pointer to inherit from
\param inputPath The path on which to start the navigation
*/
inline QString getExistingDirectory(QWidget *parent, const QString inputPath)
{
  QFileDialog fileDialog;
  fileDialog.setWindowFlags(fileDialog.windowFlags() & ~Qt::WindowContextHelpButtonHint);
  QString directory = fileDialog.getExistingDirectory(parent, "Open Directory", inputPath, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);

  if (!directory.isNull())
  {
    return directory;
  }
  else
  {
    ShowErrorMessage("No folder selected; please try again");
    return "";
  }
}

/**
\brief Open existing directory (this is the mechanism to use for saving directories as well)

\param parent The QWidget pointer to inherit from
\param inputPath The path on which to start the navigation
*/
inline QString getExistingDirectory(QWidget *parent, const std::string inputPath)
{
  return getExistingDirectory(parent, QString(inputPath.c_str()));
}

/**
\brief Open existing file

\param parent The QWidget pointer to inherit from
\param inputPath The path on which to start the navigation
\param extensions The File type and corresponding extensions as shown in the exemplary default case 
*/
inline QString getExistingFile(QWidget *parent, const QString inputPath, const QString extensions = "Images (*.nii.gz *.nii)")
{
  QFileDialog fileDialog;
  fileDialog.setWindowFlags(fileDialog.windowFlags() & ~Qt::WindowContextHelpButtonHint);
  QString filename = fileDialog.getOpenFileName(parent, "Select File", inputPath, extensions, 0, QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);

  if (!filename.isNull())
  {
    return filename;
  }
  else
  {
    ShowErrorMessage("No file selected; please try again");
    return "";
  }
}

/**
\brief Open existing file

\param parent The QWidget pointer to inherit from
\param inputPath The path on which to start the navigation
\param extensions The File type and corresponding extensions as shown in the exemplary default case
*/
inline QString getExistingFile(QWidget *parent, const std::string inputPath, const std::string extensions = "Images (*.nii.gz *.nii)")
{
  return getExistingFile(parent, QString(inputPath.c_str()), QString(extensions.c_str()));
}

/**
\brief Get file name to save (doesn't have to be existing)

\param parent The QWidget pointer to inherit from
\param inputPath The path on which to start the navigation
\param extensions The File type and corresponding extensions as shown in the exemplary default case
\param defaultFileName The file name that is populated by default
*/
inline QString getSaveFile(QWidget *parent, const QString inputPath, const QString defaultFileName = "", const QString extensions = "Images (*.nii.gz *.nii)")
{
  QFileDialog fileDialog(parent, "Save File", inputPath, extensions);
  fileDialog.setWindowFlags(fileDialog.windowFlags() & ~Qt::WindowContextHelpButtonHint);
  fileDialog.setOptions(QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);
  fileDialog.selectFile(defaultFileName);
  fileDialog.setFileMode(QFileDialog::AnyFile);
  fileDialog.setAcceptMode(QFileDialog::AcceptSave);

  //QString filename = fileDialog.getSaveFileName(parent, "Save File", inputPath.c_str(), extensions.c_str(), 0, QFileDialog::DontResolveSymlinks | QFileDialog::DontUseNativeDialog);
  int ret = fileDialog.exec();
  if (ret == QDialog::Accepted)
  {
    return fileDialog.selectedFiles()[0];
  }
  else
  {
    ShowErrorMessage("No file selected; please try again");
    return "";
  }
}

/**
\brief Get file name to save (doesn't have to be existing)

\param parent The QWidget pointer to inherit from
\param inputPath The path on which to start the navigation
\param extensions The File type and corresponding extensions as shown in the exemplary default case
\param defaultFileName The file name that is populated by default
*/
inline QString getSaveFile(QWidget *parent, const std::string inputPath, const std::string defaultFileName = "", const std::string extensions = "Images (*.nii.gz *.nii)")
{
  return getSaveFile(parent, QString(inputPath.c_str()), QString(defaultFileName.c_str()), QString(extensions.c_str()));
}
