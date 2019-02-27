# Alternative modes

This section covers methods for segmentation different than the default one. Please note that the default mode was chosen for producing the best results. Despite that, [geotrainfull](#geotrainfull) mode can also produce good results but is useful only for 2-class classification and requires thresholding to produce labels. Also, [agd](#agd) mode [(related paper)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4395536/) can produce some fast segmentations of a single region. The rest of the modes are provided mostly for documentation purposes.

### Contents 
- [Methods using just SVM](#methods-using-just-svm)
- [Methods using just AGD](#methods-using-just-agd)
- [Methods using SVM and AGD](#methods-using-svm-and-agd)
- [Methods using random forests](#methods-using-random-forests)

To specify a mode, before using ```gts.Execute();``` use 
```cpp
void SetMode(GeodesicTrainingSegmentation::MODE);
``` 

For example for mode SVM_LABELS, do:

```cpp
gts.SetMode(GeodesicTrainingSegmentation::MODE::SVM_LABELS);
```

<a name="#methods-using-just-svm"></a>
## Methods using just SVM

###### SVM_LABELS

Uses an ensemble of SVMs to segment an image. Apart from ```SetMode```, usage is the same as the default mode. In general, as with all the methods not using some sort of AGD map which provides locality to the algorithm, results tend to contain a lot of false positives.

###### SVM_PSEUDO

Only for 2 classes. There should be one label for the voxels inside the area of interest and one label for the voxels outside of it. Use whichever 2 non-zeros labels you want and provide all the images through either ```void SetInputImages(std::vector<LabelsImagePointer>);``` or ```void SetInputImages(std::vector<std::string>);```.

The output (if saving to file is enabled) is two images, one for each class, containing high values at the voxels with a high probability for being part of that class and low values for the opposite. The values range is \[0,1\]. 

Programmatically, the results can be obtained using (after ```auto executionResult = gts.Execute();```):

```cpp
typename itk::Image<float, 3>::Pointer PseudoProbImagePointer;
PseudoProbImagePointer posPseudoMap = executionResult->posImage;
PseudoProbImagePointer negPseudoMap = executionResult->negImage;
int posLabel = executionResult->posLabel; // For which label is the posPseudoMap
int negLabel = executionResult->negLabel; // For which label is the negPseudoMap
```

This method is not recommended and converting this pseudoprobability images to labels is not supported.

<a name="#methods-using-just-agd"></a>
## Methods using just AGD

<a name="#agd"></a>
###### agd

This method takes *only one* image and the "sample of labels" image as input. Samples should exist only for one label (or if there are more and you don't want to transform it, pass through ```void SetLabelOfInterest(int);``` which label to use for AGD). In other words, pass samples only for the area of interest and not for anything else. Use ```void SetThreshold(double);``` to set a custom threshold (Default is 25).

The method (if saving to file is enabled) outputs a distance map (named agd_X_classY.nii.gz). This map is then thresholded and the output is named labels_thres_classY_tZ.nii.gz.

Programmatically, the results can be obtained using (after ```auto executionResult = gts.Execute();```):

```cpp
typename itk::Image<int, 3>::Pointer AgdImagePointer;
typename itk::Image<int, 3>::Pointer LabelsImagePointer;
AgdImagePointer agdMapImage = executionResult->agdMapImage;
LabelsImagePointer labelsImage = executionResult->labelsImage;
```

A lot of the times, the threshold value used might not be the optimal one (because it is case dependent). The user can look at the output distance map image and decide that a different value should be used instead. In that case there is no need to run the algorithm again. Instead you can use:

```cpp
#include "GeodesicTrainingSegmentation.h"

...

GeodesicTrainingSegmentation::Coordinator<> gts;

std::string fullPathToDistanceMap = "C:/GeodesicTraining/agd_X_classY.nii.gz";
gts.SetInputImage(fullPathToDistanceMap); // You can also pass an itk image pointer

int newThreshold = 30; // The better threshold that was decided
gts.SetThreshold(newThreshold);

int labelToUse = 2; // Which label to put at the voxels lower than the threshold (0 is used elsewhere)
gts.SetLabelOfInterest(labelToUse);

// Optionally you can save the result to a file using SetOutputPath(std::string), SetSaveAll(bool)

// Executing 
auto executeResult = gts.Execute();

if (executeResult->ok) {
   // You can ignore labelsRenamedPtr if you only want the results saved to file
   typedef typename itk::Image< int, 3 >::Pointer   LabelsImagePointer;
   LabelsImagePointer newLabelsImage = executeResult->labelsImage; // The result segmantation mask
}
else {
   std::cout << executeResult->errorMessage << "\n";
}
```

<a name="#methods-using-svm-and-agd"></a>
## Methods using SVM and AGD

###### geotrain

This method is not recommended. For distance maps it is recommended to use [geotrainfull](#geotrainfull) for better results or [agd](#agd) for faster results.

Only for 2 classes. There should be one label for the voxels inside the area of interest and one label for the voxels outside of it. Provide all the images through  either ```void SetInputImages(std::vector<LabelsImagePointer>);``` or ```void SetInputImages(std::vector<std::string>);``` and use whichever 2 non-zeros labels you want. Specify which of the two labels is the label for the area of interest using ```void SetLabelOfInterest(int);```. Default label of interest is 1.

The methods (if saving to file is enabled) outputs two distance maps (named agd_X_classY.nii.gz). The map for the LABEL_OF_INTEREST is then thresholded with the value provided in ```void SetThreshold(double);``` (the output labels image is named labels_thres_classY_tZ.nii.gz). Default threshold value if the parameter is not provided is 25, but 25 is not recommended for this mode, a value like 100 is better. 

Programmatically, the results can be obtained using (after ```auto executionResult = gts.Execute();```):

```cpp
typename itk::Image<int, 3>::Pointer AgdImagePointer;
typename itk::Image<int, 3>::Pointer LabelsImagePointer;
AgdImagePointer agdMapImage = executionResult->agdMapImage;
LabelsImagePointer labelsImage = executionResult->labelsImage;
```

A lot of the times, the threshold value used might not be the optimal one (because it is case dependent). The user can look at the output distance map image and decide that a different value should be used instead. In that case there is no need to run the algorithm again. Instead you can use:

```cpp
#include "GeodesicTrainingSegmentation.h"

...

GeodesicTrainingSegmentation::Coordinator<> gts;

std::string fullPathToDistanceMap = "C:/GeodesicTraining/agd_X_classY.nii.gz";
gts.SetInputImage(fullPathToDistanceMap); // You can also pass an itk image pointer

int newThreshold = 30; // The better threshold that was decided
gts.SetThreshold(newThreshold);

int labelToUse = 2; // Which label to put at the voxels lower than the threshold (0 is used elsewhere)
gts.SetLabelOfInterest(labelToUse);

// Optionally you can save the result to a file using SetOutputPath(std::string), SetSaveAll(bool)

// Executing 
auto executeResult = gts.Execute();

if (executeResult->ok) {
   // You can ignore labelsRenamedPtr if you only want the results saved to file
   typedef typename itk::Image< int, 3 >::Pointer   LabelsImagePointer;
   LabelsImagePointer newLabelsImage = executeResult->labelsImage; // The result segmantation mask
}
else {
   std::cout << executeResult->errorMessage << "\n";
}
```


where LABEL_OF_INTEREST is which value to put where the distance map has values lower than the threshold (0 is used elsewhere)

<a name="#geotrainfull"></a>
###### geotrainfull

Only for 2 classes. There should be one label for the voxels inside the area of interest and one label for the voxels outside of it. Provide all the images through  either ```void SetInputImages(std::vector<LabelsImagePointer>);``` or ```void SetInputImages(std::vector<std::string>);``` and use whichever 2 non-zeros labels you want. Specify which of the two labels is the label for the area of interest using ```void SetLabelOfInterest(int);```. Default label of interest is 1.

The methods (if saving to file is enabled) outputs two distance maps (named agd_X_classY.nii.gz). The map for the LABEL_OF_INTEREST is then thresholded with the value provided in ```void SetThreshold(double);``` (the output labels image is named labels_thres_classY_tZ.nii.gz). Default threshold value if the parameter is not provided is 25, but 25 is not recommended for this mode, a value like 100 is better. 

Programmatically, the results can be obtained using (after ```auto executionResult = gts.Execute();```):

```cpp
typename itk::Image<int, 3>::Pointer AgdImagePointer;
typename itk::Image<int, 3>::Pointer LabelsImagePointer;
AgdImagePointer agdMapImage = executionResult->agdMapImage;
LabelsImagePointer labelsImage = executionResult->labelsImage;
```

A lot of the times, the threshold value used might not be the optimal one (because it is case dependent). The user can look at the output distance map image and decide that a different value should be used instead. In that case there is no need to run the algorithm again. Instead you can use:

```cpp
#include "GeodesicTrainingSegmentation.h"

...

GeodesicTrainingSegmentation::Coordinator<> gts;

std::string fullPathToDistanceMap = "C:/GeodesicTraining/agd_X_classY.nii.gz";
gts.SetInputImage(fullPathToDistanceMap); // You can also pass an itk image pointer

int newThreshold = 30; // The better threshold that was decided
gts.SetThreshold(newThreshold);

int labelToUse = 2; // Which label to put at the voxels lower than the threshold (0 is used elsewhere)
gts.SetLabelOfInterest(labelToUse); 

// Optionally you can save the result to a file using SetOutputPath(std::string), SetSaveAll(bool)

// Executing 
auto executeResult = gts.Execute();

if (executeResult->ok) {
   // You can ignore labelsRenamedPtr if you only want the results saved to file
   typedef typename itk::Image< int, 3 >::Pointer   LabelsImagePointer;
   LabelsImagePointer newLabelsImage = executeResult->labelsImage; // The result segmantation mask
}
else {
   std::cout << executeResult->errorMessage << "\n";
}
```

where LABEL_OF_INTEREST is which value to put where the distance map has values lower than the threshold (0 is used elsewhere)

<a name="#methods-using-random-forests"></a>
## Methods using random forests

Generally Random Forests where found to be worse than SVMs for the task.

<a name="#rf"></a>
###### rf

The simplest random forest mode.

Usage is the same as the default mode, except that you can optionally pass custom parameters for the random forest using ```void SetRfConfigFile(std::string);```. [See here](../extra/configurations/rf/rf_example_conf.config) for an example custom RF configuration.

In general, as with all the methods not using some sort of AGD map which provides locality to the algorithm, results tend to contain a lot of false positives.

###### rfauto

Usage is the same as the default mode. Difference with [rf](#rf) is that the parameters are tuned automatically.

In general, as with all the methods not using some sort of AGD map which provides locality to the algorithm, results tend to contain a lot of false positives.

<a name="#agdrf"></a>
###### agdrf

Uses AGD maps alongside the (unprocessed) input images as input to the RF.

Usage is the same as the default mode, except that you can optionally pass custom parameters for the random forest using ```void SetRfConfigFile(std::string);```. [See here](../extra/configurations/rf/rf_example_conf.config) for an example custom RF configuration.

Uses AGD maps alongside the (unprocessed) input images as input to the RF. 

###### agdrfauto

Uses AGD maps alongside the (unprocessed) input images as input to the RF.

Usage is the same as the default mode. Difference with [agdrf](#agdrf) is that the parameters are tuned automatically.
