# Advanced Usage

### Contents 
- [Use different labels](#use-different-labels)
- [Change labels in a labels image](#change-labels-in-a-labels-image)
- [Compare segmentation to ground truth](#compare-segmentation-to-ground-truth)
- [More optional parameters](#more-optional-parameters)
- [Isolate an area in an image](#isolate-an-area-in-an-image)
- [Custom SVM ensemble](#custom-svm-ensemble)

<a name="#use-different-labels"></a>
## Use different labels

If you want to use different labels than the ones specified [here](README.md#input-labels) use methods:
```cpp
void SetTumorCoreLabelMRI(int labelTC);
```
```cpp
void SetEnhancedTumorLabelMRI(int labelET);
```
```cpp
void SetEdemaLabelMRI(int labelED);
```
```cpp
void SetHealthyTissueLabelMRI(int labelHT);
``` 

where labelTC/labelET/labelED/labelHT are non-zero integers.

<a name="#change-labels-in-a-labels-image"></a>
## Change labels in a labels image

```cpp
#include "GeodesicTrainingSegmentation.h"
...
typedef typename itk::Image< int, 3 >::Pointer LabelsImagePointer;
LabelsImagePointer labelsPtr = ...; 

std::unordered_map<int, int> changeLabelsMap;

changeLabelsMap[2] = 1; // Will change label 2 to 1
changeLabelsMap[3] = 0; // Will change label 3 to 0

GeodesicTrainingSegmentation::Coordinator<> gts();
gts.SetMode(GeodesicTrainingSegmentation::MODE::CHANGE_LABELS);
gts.SetLabels(labelsPtr);
gts.SetChangeLabelsMap(changeLabelsMap);

auto executeResult = gts.Execute();

if (executeResult->ok) {
   LabelsImagePointer labelsRenamedPtr = executeResult->labelsImage;
}
else {
   std::cout << executeResult->errorMessage << "\n";
}

```

```SetChangeLabelsMap``` can also be used in all the modes, like the code described [here](README.md/#run). Note if ```SetChangeLabelsMap``` is not used, labels for healthy tissue will be removed. Pass an empty changeLabelsMap, if you don't want that. 

<a name="#compare-segmentation-to-ground-truth"></a>
## Compare segmentation to ground truth

```cpp
#include "GeodesicTrainingSegmentation.h"
...
typedef typename itk::Image< int, 3 >::Pointer LabelsImagePointer;
LabelsImagePointer labelsPtr = ...; 
LabelsImagePointer groundTruthPtr = ...;

GeodesicTrainingSegmentation::Coordinator<> gts();
gts.SetMode(GeodesicTrainingSegmentation::MODE::CHECK_ACCURACY);
gts.SetLabels(labelsPtr); // The segmentation to be compared to ground truth
gts.SetGroundTruth(groundTruthPtr);

auto executeResult = gts.Execute();


if (executeResult->ok) {
   double diceScore = executeResult->diceScoreAll;
   double sensitivity = executeResult->sensitivityAll;
   int falsePositives = executeResult->falsePositiveCountAll;
}
else {
   std::cout << executeResult->errorMessage << "\n";
}
```

If you want some labels to be skipped when comparing with a ground truth segmentation, change before calling ```Execute```:

```cpp
std::vector<int> groundTruthSkipLabels;
groundTruthSkipLabels.push_back(1); // Will skip label 1 when comparing

gts.SetGroundTruth(groundTruthPtr, groundTruthSkipLabels);

```

Please note that label 0 is the null label and is not used at the calculations regardless.

Additionally, the ```executeResult``` contains dice scores, sensitivity and false positives for each individual label used.

```cpp
std::map< int, double> diceScoreForLabel = executeResult->diceScore;
std::map< int, double> sensitivityForLabel = executeResult->sensitivity;
std::map< int, int>    falsePositives = executeResult->falsePositivesCount;

// For example the diceScore for label 2 can be found doing:
double diceScoreForLabel2 = diceScoreForLabel[2];
```

Comparing to a ground truth segmentation can also be used in all the modes, like the code described [here](README.md/#run).

<a name="#more-optional-parameters"></a>
## More optional parameters

* ```SetTimerEnabled(bool)``` Records duration and saves it to "time_report.txt"
* ```SetNumberOfThreads(int)``` Sets the number of threads manually. Default is 16.

To also save the results to a file use both parameters:

```cpp
gts.SetOutputPath(fullPathToOutputFolderAsStdString);
gts.SetSaveAll(true);
```

<a name="#isolate-an-area-in-an-image"></a>
## Isolate an area in an image

```cpp
#include "GeodesicTrainingSegmentation.h"
...
typedef typename itk::Image< float, 3 >::Pointer InputImagePointer;
typedef typename itk::Image< int, 3 >::Pointer LabelsImagePointer;
LabelsImagePointer labelsPtr = ...; // The labels here will be used
InputImagePointer inputImagePtr = ...; // The MRI image

// This will keep only the voxels from inputImagePtr
// that have the label 2 in labelsPtr
int labelOfInterest = 2;

GeodesicTrainingSegmentation::Coordinator<> gts();
gts.SetMode(GeodesicTrainingSegmentation::MODE::SEGMENT);
gts.SetInputImage(inputImagePtr);
gts.SetLabels(labelsPtr);
gts.SetLabelOfInterest(labelOfInterest);

auto executeResult = gts.Execute();

if (executeResult->ok) {
   InputImagePointer isolatedImage = executeResult->segmentedFloatImage;
}
else {
   std::cout << executeResult->errorMessage << "\n";
}
```

The output will be the input image where the voxels that *don't* have the label specified in variable labelOfInterest become zero.

<a name="#custom-svm-ensemble"></a>
## Custom SVM ensemble

By default no configuration file need to be provided and an ensemble of SVMs with RBF, chi2 and histogram intersection kernels are trained, but a custom YAML configuration file can be provided using ```void SetConfigFile(std::string);```

<a name="#the-different-keywords-are"></a>
##### The different keywords are:

* __svms__ The root node should always be "svms"
* __kernel_type__ Different kernels available are: linear, rbf, poly, sigmoid, chi2, inter
* __type__ If unsure of what SVM type is remove it, or leave it as C_SVC. See OpenCV's docs for more info.
* __kfold__ If unsure of what kfold is remove it, or leave it at 10 (default) or 5 (faster)
* __consider_weights__ Misclassifying a class that has less label samples is penalized more during training. C_SVC only.
* __importance__ Remove it for using importance=1 (default)
  * In case of normal classification, each vote from a SVM is weighted by its importance value. For example in case of 3 SVMs if SVM_A has importance 0.4 and SVM_B, SVM_C both have 0.3 then the vote of SVM_A is always used as output, except if SVM_B and SVM_C both vote for the same different class.
  * In the modes that produce pseudoprobability maps the average of all the pseudoprobabilities of the different SVMs is taken but it can also be weighted using importance so that higher importance impacts the result more.
* __term_criteria__ Remove it for default term criteria
  * __criteria_type__ Can be "MAX_ITER", or "EPS" or "MAX_ITER+EPS"
  * __max_count__ Max iterations
  * __epsilon__ Max error
* __pretrained__ Use a pretrained OpenCV svm model that was saved to a file
* __c/gamma/p/nu/coef/degree__ Values/Ranges for different parameters (See OpenCV docs for more info). A default value is used if one of this keywords is missing. Different options are:
  * *Specific value*
  * __auto__ For using the default grid of the parameter for optimization.
  * __min_value,max_value,log_step__ {min_value, min_value * log_step, min_value * log_step^2, ...} will be tried, while min_value * log_step ^ n < max_value.  

The different SVM configurations start with "-" (list items). Example configuration files can be found at the [extra/configurations/svm](../extra/configurations/svm) directory. The following configuration is for demonstration purposes and not to be used as input.

```yaml
---
# This is a comment
svms:
  - kernel_type: linear
    kfold: 5
    consider_weights: false 
    c: 0.8
  - kernel_type: rbf
    type: C_SVC
    kfold: 10 
    importance: 0.4
    consider_weights: true 
    c: 0.7
    gamma:
      min_value: 0.001 
      max_value: 1.0 
      log_step: 2.0 
  - kernel_type: sigmoid
    consider_weights: true 
    term_criteria:
      criteria_type: MAX_ITER
      max_count: 900
      epsilon: 0.01
    c: 0.9
    gamma: auto
    coef: auto
  - pretrained: C:\GeodesicTraining\model_linear1.xml
    importance: 0.2
...
``` 
