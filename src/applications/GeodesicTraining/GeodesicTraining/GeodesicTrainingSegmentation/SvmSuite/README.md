# SvmSuite - A configurable library for SVM operations

Depends on OpenCV and uses OpenCV's matrices (cv::Mat) for input/output. Requires C++11.

## Build

OpenCV needs to be installed in your system and the path to OpenCV needs to be provided to cmake. Assuming the directory SvmSuite is copied inside your project's root directory, an example CMakeLists.txt file could look like this:

```cmake
cmake_minimum_required(VERSION 3.0)
SET( PROJECT_NAME YourProjectName)
SET( PROJECT_VERSION "1.0")
ADD_DEFINITIONS(-DPROJECT_VERSION="${PROJECT_VERSION}" )

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

# For OpenCV
FIND_PACKAGE( OpenCV REQUIRED )

# Setting your sources
FILE(GLOB_RECURSE Sources "${PROJECT_SOURCE_DIR}/*.*")

# Including the paths for headers
include_directories(
  "${PROJECT_SOURCE_DIR}/SvmSuite/include/SvmSuite"
)

# Adding the subproject
add_subdirectory(SvmSuite)

# Add sources to executable
ADD_EXECUTABLE( ${PROJECT_NAME} 
  ${sources}
)

# The paths to the API of the subproject
target_include_directories( ${PROJECT_NAME}
  PUBLIC SvmSuite/include/SvmSuite
)

# Link the libraries to be used
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
  SvmSuite
  ${OpenCV_LIBRARIES}
)
```

## Usage

Execution happens through the ```SvmSuite::Manager``` class.

Also, see file ```SvmSuiteDescription.h``` for details about ```SvmSuite::SvmDescription```. In short, a single ```SvmSuite::SvmDescription``` represents a configuration for a single SVM model. You can use an ensemble of different SVM models for making predictions. That is achieved by passing a vector of SvmDescription(s) to ```SvmSuite::Manager```. If you want to use the default configuration (an ensemble of rbf, chi2 and histogram intersection kernels where all parameters are tuned automatically for all three kernels) get it through:

```cpp
std::vector<SvmSuite::SvmDescription> svmDescs = SvmSuite::SvmDescription::GetDefaultSvmDescriptions();
```

#### Example

```cpp
#include "SvmSuite.h"

...

// The training data should be float (CV_32F type). The labels data should be int (CV_32S type).

cv::Mat trainingData = ...; // X rows, Y cols. One row per sample. Each column is a feature.
cv::Mat labelsMat = ...; // X rows, 1 cols. The responses for each row of training data

// Optionally provide a weights mat.
// Let Z be the number of unique different responses (labels).
// Weights mat will should have Z rows, 1 cols.
// Imagine the labels are sorted in ascending order. (according to their value, not their sample size)
// The first weight should be for the smaller label, the last for the bigger one.
// For example suppose you use labels: 2, -1, 4. Row 1 should have the weight of label -1.
cv::Mat weightsMat = ...; 

SvmSuite::Manager svmManager;

std::vector<SvmSuite::SvmDescription> svmDescs = SvmSuite::SvmDescription::GetDefaultSvmDescriptions();
svmManager.AddSvmDescriptions(svmDescs);

svmManager.SetTrainData(trainingMat, labelsMat);

// If you want to pass weights for each class use instead: 
// svmManager.SetTrainData(trainingMat, labelsMat, weightsMat);

// More Optional Methods
svmManager.SetOutputPath("C:/GeodesicTraining"); // Folder to optionally write results
svmManager.SetTimerEnabled(true); // Will create time_report.txt in output path
svmManager.SetSavingModelsEnabled(true); // Will save the trained model(s) in output path
svmManager.SetVerbose(true); // Show messages in the console 
svmManager.SetNumberOfThreads(32); // 32 is the default

svmManager.Train();

bool pseudoProbMapResult = false;

// Change to true for probabilities (instead of responses)
// Probabilities is 2-class only. Responses is n-class.

// Testing mat contains the samples to test
// W rows, Y cols. One row per sample. Each column is a feature.
// There should exactly the same number of features for trainingMat and testingMat
cv::Mat testingMat = ...; 

auto svmTestingResult = svmManager.Test(testingMat, predictFlag);

if (predictFlag) {
  // Probabilities
  // The posMat contains probabilities for each sample to belong to class posLabel
  // The reverse for negMat (probabilities are when using 2-class classification only)
  cv::Mat posMat = svmTestingResult->posMat;
  cv::Mat negMat = svmTestingResult->negMat;
  int posLabel = svmTestingResult->posLabel; // Which of the two labels is posMat about
  int negLabel = svmTestingResult->negLabel;
}
else {
  // Responses
  cv::Mat responsesResult = svmTestingResult->labelsMat; // One response per row
}
```

## Documentation

SvmSuite::Manager Public Methods:

- Adders and Setters

```cpp
void AddSvmDescriptions(std::vector< SvmSuite::SvmDescription > svmDescs);
void AddPretrainedModel(std::string pretrainedModelPath);
void AddSvmsFromConfig(std::string configPath);
void SetTrainData(cv::Mat &trainingMat, cv::Mat &labelsMat, cv::Mat &weightsMat);
void SetVerbose(bool verbose);
void SetOutputPath(std::string path);
void SetSavingModelsEnabled(bool modelsEnabled);
void SetTimerEnabled(bool timerEnabled);
void SetNumberOfThreads(int numberOfThreads);
```

- The result returned after testing the model(s) is a ```std::shared_ptr``` of:

```cpp
typedef struct Result {
  cv::Mat posMat;
  cv::Mat negMat;
  int posLabel;
  int negLabel;
  cv::Mat labelsMat;
} Result;
```

- Training and Testing

```cpp
/**
Train the models specified in the svm descriptions
*/
void Train();

/**
Test using the ensemble of trained svm models
@param testingMat cv::Mat where the rows are the samples to test and columns are the features
@param pseudoProbMapResult true->  the result would be pos and neg pseudoprobability maps (only for n=2 classification)
                           false-> the result would be the predicted labels (n>=2 classification)
@return pointer to SvmSuite::Manager::Result object (which contains the result images)
*/
std::shared_ptr<Result> Test(cv::Mat &testingMat, bool pseudoProbMapResult = false);

/**
Test using the ensemble of trained svm models
@param testingMat cv::Mat where the rows are the samples to test and columns are the features
@param skipZerosMat cv::Mat that has as many rows as testingMat, 
                    where for each row if skipZerosMat has zero the line will be skipped
                    and the label 0 will be set. 
                    In practice you can pass testingMat two times (extra columns would not matter)
@param pseudoProbMapResult true->  the result would be pos and neg pseudoprobability maps (only for n=2 classification)
                           false-> the result would be the predicted labels (n>=2 classification)
@param skipZeros whether to user skipZerosMat
@return pointer to SvmSuite::Manager::Result object (which contains the result images)
*/
std::shared_ptr<Result> Test(cv::Mat &testingMat, cv::Mat &skipZerosMat, bool pseudoProbMapResult = false, bool skipZeros = true);
```
