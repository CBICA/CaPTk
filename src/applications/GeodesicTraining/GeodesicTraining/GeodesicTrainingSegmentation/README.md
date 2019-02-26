# Geodesic Training - Tumor Segmentation (Library)

This library is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania.

### Contents 
- [What can it do?](#what-can-it-do)
- [How to build?](#how-to-build)
- [How to use?](#how-to-use)


<a name="what-can-it-do"></a>
## What can it do?

Produce masks for the different areas specified in the "sample of labels" image.

For brain tumors specifically, produce masks for:
- necrotic and non-enhancing tumor core
- GD-enhancing tumor 
- peritumoral edema
- other areas

See [advanced usage](EXTRA_ADVANCED.md) section for complementary operations.

See [alternative modes](EXTRA_MODES.md) section for alternative methods for segmentation.

<a name="how-to-build"></a>
## How to build

Depends on ITK, OpenCV. Requires a C++11 compliant compiler and cmake.

Assuming the GeodesicTrainingSegmentation directory has been copied inside your project's root directory and your sources are at src directory, an example CMakeLists.txt file at your the root directory of your project could look like this:

```cmake
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

SET( PROJECT_NAME YourProjectName )
PROJECT( ${PROJECT_NAME} )
 
# For ITK
FIND_PACKAGE( ITK REQUIRED )
SET(ITK_NO_IO_FACTORY_REGISTER_MANAGER "OFF")
INCLUDE( ${ITK_USE_FILE} )

# For OpenCV
FIND_PACKAGE( OpenCV REQUIRED )
#INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

# For C++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 

# For OpenMP
FIND_PACKAGE(OpenMP REQUIRED)
SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

# Finding your source files
FILE( GLOB_RECURSE sources "${PROJECT_SOURCE_DIR}/src/*.*" )

# Including the paths for headers
INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/GeodesicTrainingSegmentation/include/GeodesicTrainingSegmentation
)

# Adding the subproject
add_subdirectory(GeodesicTrainingSegmentation)

# Add your sources to the executable
ADD_EXECUTABLE( ${PROJECT_NAME} 
  ${sources}
)

# Link the libraries to be used
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
  ${ITK_LIBRARIES}
  ${OpenCV_LIBRARIES}
  GeodesicTrainingSegmentation
)
```

ITK and OpenCV should be installed on your system. You have to pass the paths to ITK and OpenCV include directories to cmake. 

###### Linux

An example cmake command on linux could look like this (but the CMAKE_INSTALL_PREFIX and the ITK version might be different):

```console
cmake -DCMAKE_INSTALL_PREFIX="/usr/local" -DCMAKE_PREFIX_PATH="${CMAKE_INSTALL_PREFIX}/include/ITK-4.13;${CMAKE_INSTALL_PREFIX}/include/opencv4" -DCMAKE_MODULE_PATH="${CMAKE_INSTALL_PREFIX}/include/opencv4" ..
```

For more information on installing the dependencies on linux you can see [here for ITK](https://itk.org/Wiki/ITK_Configuring_and_Building_for_Ubuntu_Linux) and [here for OpenCV](https://docs.opencv.org/3.4/d7/d9f/tutorial_linux_install.html), but you can also try to look for newer versions.

###### Windows

This library can be built on Windows by providing the paths to ITK and OpenCV. See [here for ITK](https://itk.org/ITK/resources/software.html) and [here for OpenCV](https://docs.opencv.org/2.4/doc/tutorials/introduction/windows_install/windows_install.html), but you can also try to look for newer versions.

In general, the ```CMAKE_PREFIX_PATH``` should look something like this ```"C:/Libraries/InsightToolkit-4.13.1/build;C:/Libraries/opencv/build/x64/vc15/lib"``` and ```CMAKE_MODULE_PATH``` something like this ```"C:/Libraries/InsightToolkit-4.13.1/build/lib/cmake/ITK-4.13"```

If when trying to run your application you get prompted that ```opencv_world343.dll``` or  ```opencv_world343d.dll``` is missing, manually copy it from opencv's build path to the directory where the your executable is located. (The path to the dll will look something like this ```C:\Libraries\opencv\build\x64\vc15\bin```)

## How to use?

#### Input files

Supports 3D NIfTI and DICOM images (.nii, .nii.gz, .dcm, .dicom files). Input files can be set either with a std::string to the full path of the image, or by a pointer to an ```itk::Image``` (see [run section](#run) below for an example). Only 3D images are supported for now.

###### Input Images

The images can be MR, CT, PET or anything else (or a combination of all). The restriction is that all images should have exactly the same size and be co-registered images of one subject taken during a single session. 

Either set by using a std::string to the full path of the image or by a pointer to an ```itk::Image< float, 3 >```. (If your input image is of another type like ``` itk::Image< int, 3 >``` you can use ITK's CastImageFilter to convert it to float, or pass it using its path).

Special operations are performed for the modalities FLAIR, T1, T1CE (Gd) and T2 when segmenting brain tumors. Although, none of them are mandatory and other modalities can be provided too.

<a name="#input-labels"></a>
###### Input Labels

The sample of labels image should also be an image that has exactly the same size as the different input images. 

The unlabeled voxels should have the value 0. Any other non-zero integer label can be used for the different areas. There should be at least two different labels.

<a name="#brain-tumor"></a>
###### Brain tumor

In the special case of brain tumors, the values of the voxels of this image should be:
* __0__  At the unlabeled voxels
* __1__  At the necrotic and non-enhancing tumor core samples
* __2__  At the peritumoral edema samples
* __3__  At the healthy tissue samples
* __4__  At the GD-enhancing tumor samples
* __5+__ At other areas the user wants to mark

Please note that *none* of the brain tumor specific classes (1,2,3,4) are mandatory even when segmenting brain tumors, but there should be at least 2 different classes in the image (there could all be 5+). Also note that apart from the expect output mask, a second one will be created where the label for healthy tissue is removed.

Either set by using a std::string to the full path of the image or by a pointer to an ```itk::Image< int, 3 >```. Please notice that the input images should have __float__ pixel type, while the labels image should have __int__. (If your input labels are of another type like ``` itk::Image< short, 3 >``` you can use ITK's CastImageFilter to convert it to int, or pass it using its path).

See [here](EXTRA_ADVANCED.md#use-different-labels) for instructions on how to use different labels for the tumor-specific classes. 

Please also note that the quality of the segmentations are highly correlated to the quality of the provided "sample of labels" image. Also, especially in the more difficult cases, there might be misclassified regions in the output segmentation. If that is the case the user can add more labels in the "sample of labels" images that are within these regions and run the tool again. 

[CaPTk](https://www.med.upenn.edu/cbica/captk/) can be used to draw the labels image. Drag and drop the MRI images in the "Images" area, draw labels through the "Drawing" tab and export using File>Save>ROI. For drawing more labels to an old labels image use File>Load>ROI, draw the new labels and export.

<a name="#run"></a>
#### Run

All the interaction is done using the class
```
GeodesicTrainingSegmentation::Coordinator<typename PixelType = float, unsinged int Dimensions = 3>
```

Please note that the template refers to the input images. Right now only float PixelType is supported, but it is kept in case there is reason to support more in the future. Dimensions can be 2 or 3. The labels image should be itk::Image<int, Dimensions>.

*Note:* If you are segmenting brain tumors and at least one of the FLAIR/T1/T1CE/T2 MRI modalities are present, pass these modalities using one of the:

```
// For 2D
void GeodesicTrainingSegmentation::Coordinator<float,2>::SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI, std::string);
void GeodesicTrainingSegmentation::Coordinator<float,2>::SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI, typename itk::Image<float, 2>::Pointer);

// For 3D
void GeodesicTrainingSegmentation::Coordinator<float,3>::SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI, std::string);
void GeodesicTrainingSegmentation::Coordinator<float,3>::SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI, typename itk::Image<float, 3>::Pointer);
```

Pass other images (either in the general non brain tumor case, or extra images that are not FLAIR/T1/T1CE/T2 for brain tumors using:

```
// For 2D
void GeodesicTrainingSegmentation::Coordinator<float,2>::SetInputImages(GeodesicTrainingSegmentation::MODALITY_MRI, std::vector<std::string>);
void GeodesicTrainingSegmentation::Coordinator<float,2>::SetInputImages(GeodesicTrainingSegmentation::MODALITY_MRI, std::vector< typename itk::Image<float, 2>::Pointer >);

// For 3D
void GeodesicTrainingSegmentation::Coordinator<float,3>::SetInputImages(GeodesicTrainingSegmentation::MODALITY_MRI, std::vector<std::string>);
void GeodesicTrainingSegmentation::Coordinator<float,3>::SetInputImages(GeodesicTrainingSegmentation::MODALITY_MRI, std::vector< typename itk::Image<float, 3>::Pointer >);

```

You can supply some FLAIR/T1/T1CE/T2 images using paths and some using pointers. ```SetInputImages``` will not override the changes made by ```SetInputImageMRI```. 

Either supply all the other images using a vector of paths or a vector of pointers. Don't call ```SetInputImages``` twice.

In the following example all images are 3D and
- flair image is passed by it's full path
- t1ce image is passed using a pointer to an itk::Image<float, 3>
- 2 other images are passed using their full path
- the "sample of labels" image is passed using a pointer to an itk::Image<int, 3>

```cpp
#include "GeodesicTrainingSegmentation.h"
...

typedef itk::Image<float, 3> InputImageType;
typedef itk::Image<int, 3>   LabelsImageType;
typedef typename InputImageType::Pointer  InputImagePointer;
typedef typename LabelsImageType::Pointer LabelsImagePointer;

std::string flairImageFullPath  = "C:/GeodesicTraining/flair.nii.gz"; // The flair MRI image
std::string otherImage1FullPath = "C:/GeodesicTraining/other_image1.nii.gz"; // Another image
std::string otherImage2FullPath = "C:/GeodesicTraining/other_image2.nii.gz"; // Another image

std::string t1ceImageFullPath  = "C:/GeodesicTraining/t1ce.nii.gz"; // The t1ce MRI image
std::string labelsFullPath     = "C:/GeodesicTraining/mask.nii.gz"; // The sample of labels image

typedef itk::ImageFileReader<InputImageType>  InputReaderType;
typedef itk::ImageFileReader<LabelsImageType> LabelsReaderType;

typename ReaderType::Pointer inputReader  = InputReaderType::New();
typename ReaderType::Pointer labelsReader = LabelsReaderType::New();

inputReader->SetFileName(t1ceImageFullPath);
InputImagePointer t1ceImage = inputReader->GetOutput();

labelsReader->SetFileName(labelsFullPath);
LabelsImagePointer labelsImage = labelsReader->GetOutput();

std::vector<std::string> otherInputImagesPaths = { otherImage1FullPath, otherImage2FullPath }; 


// All the interaction happens through the coordinator class
GeodesicTrainingSegmentation::Coordinator<float, 3> gts;

// Setting the parameters
gts.SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI::FLAIR, flairImageFullPath);
gts.SetInputImageMRI(GeodesicTrainingSegmentation::MODALITY_MRI::T1CE,  t1ceImage);
gts.SetInputImages(otherInputImagesPaths); // Will not override FLAIR, T1, T1CE or T2 
gts.SetLabels(labelsImage);

// Also optionally you can save the result to a file using SetOutputPath(std::string), SetSaveAll(bool)
std::string outputDirWhichWillContainTheResult = "C:/GeodesicTraining/output";
gts.SetOutputPath(outputDirWhichWillContainTheResult);
gts.SetSaveAll(true);

// Executing 
auto executeResult = gts.Execute();

if (executeResult->ok) {
   // You can ignore labelsRenamedPtr if you only want the results saved to file
   LabelsImagePointer labelsRenamedPtr = executeResult->labelsImage; // The result segmantation mask
}
else {
   std::cerr << executeResult->errorMessage << "\n";
}
```

For more advanced usage (other optional parameters, using different labels, comparing segmentation to ground truth, optional parameters, isolating areas in an image and custom SVM configurations) see the [advanced usage](EXTRA_ADVANCED.md) section.
