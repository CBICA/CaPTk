#pragma once
#include "cbicaUtilities.h"
#include "itkVariableLengthVector.h"
#include "itkVariableSizeMatrix.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"


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

#define USE_LANDMARKS_RAS_COORD // [TBD] - this isn't needed since we will be working with a single coordinate system in the future

#define SVM_MODEL_EXTENSIONS "Models (*.model)"


#define NII_EXT ".nii" // [TBD] - convert to enum
#define NII_GZ_EXT ".nii.gz"
#define HDR_EXT ".hdr"
#define IMG_EXT ".img"
#define PARAM_EXT ".txt"
#define CSV_EXT ".csv"

//  static QString IMAGES_EXTENSIONS = "Images (*.nii.gz *.nii *.dcm)";

  const std::string captk_currentApplicationPath = cbica::getExecutablePath();
  const std::string loggerFolderBase = cbica::getUserHomeDirectory() + "/." + std::string(PROJECT_NAME) + "/"; // kept separate because of folder creation
  const std::string captk_StuffFolderBase = cbica::getUserHomeDirectory() + "/" + std::string(PROJECT_NAME) + "/"; // kept separate because of folder creation
  const std::string captk_SampleDataFolder = cbica::getUserHomeDirectory() + "/" + std::string(PROJECT_NAME) + "/SampleData"; // kept separate because of folder creation
  const std::string captk_PretrainedFolder = cbica::getUserHomeDirectory() + "/" + std::string(PROJECT_NAME) + "/PretrainedModels"; // kept separate because of folder creation
  const std::string loggerFolder = loggerFolderBase + std::string(PROJECT_VERSION) + "/";

  const std::string loggerFile = loggerFolder + cbica::replaceString(cbica::getCurrentLocalDate(), ":", "-") + ".log";
  const std::string aboutScreenSeen = loggerFolder + "aboutSeen.txt"; // file to check if the user has seen "about" page or not
  const std::string tutorialScreen = loggerFolder + "tutorial.txt"; // file to check if the user has seen "about" page or not
  const std::string closeConfirmation = loggerFolder + "close.txt"; // file to check if the user has seen "about" page or not

  // Common defines - these are required for GCC compilation 
  using VariableLengthVectorType = itk::VariableLengthVector< double >; // TBD: move from itk::matrix to vnl::matrix for simplicity
  using VariableSizeMatrixType = itk::VariableSizeMatrix< double >; // TBD: move from itk::matrix to vnl::matrix for simplicity
  using VectorVectorDouble = std::vector< std::vector < double > >;
  using VectorDouble = std::vector < double >;
  using ImageTypeFloat3D = itk::Image< float, 3 >;
  using ImageTypeFloat4D = itk::Image< float, 4 >;
  using ImageTypeShort3D = itk::Image< short, 3 >;
  using ImageTypeFloat3DIterator = itk::ImageRegionIteratorWithIndex< ImageTypeFloat3D >;