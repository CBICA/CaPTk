# Geodesic Training - Semantic Segmentation

<!--[![Version](https://img.shields.io/github/release/CBICA/GeodesicTraining.svg)](https://github.com/CBICA/GeodesicTraining/releases)-->

A semi-automatic, configurable tool that uses an ensemble of Support Vector Machines and Adaptive Geodesic Distance to segment regions of interest in medical images. [[Download](https://github.com/CBICA/GeodesicTraining/releases)]

This tool is developed and maintained by the [Center for Biomedical Image Computing and Analytics (CBICA)](https://www.cbica.upenn.edu/) at the University of Pennsylvania.


### Disclaimer
- The software has been designed for research purposes only and has neither been reviewed nor approved for clinical use by the Food and Drug Administration (FDA) or by any other federal/state agency.
- This code (excluding dependent libraries) is governed by the license provided [here](http://www.med.upenn.edu/sbia/software-agreement.html) unless otherwise specified.

### Contents 
- [How does it work?](#how-does-it-work)  
- [What can it do?](#what-can-it-do)
- [How to use?](#how-to-use)
- [How to build?](#how-to-build)
- [How to use as a library?](#how-to-use-as-a-library)

<a name="how-does-it-work"></a>
## How does it work?

During an MRI scan, or a similar procedure like CT or PET, multiple images of the same patient are captured using different modalities (T1, T2, etc). This tool uses these images along with a (small and quickly drawn) sample of labels to produce a segmentation of different tumor areas.

#### Algorithm

The algorithm works individually for each patient (no pretrained models are used). The input is the images taken with different modalities for the patient and a sample of labels. All the input images have exactly the same size.

First operation is the AGD algorithm (described below). For each input image one or more distance maps are produced, depending on the the modality. Second operation is training an ensemble of SVM classifiers on the spot using as input the (unprocessed) images and the distance maps. The voxel intensities of the different images and maps act as features for training. The trained model is then used to make predictions for the unlabeled regions.

###### AGD (Adaptive Geodesic Distance) algorithm

The algorithm takes as input one image and the "sample of labels" image with all the samples removed except for the ones belonging to one of the classes. The algorithm uses the geodesic distance and the intensity changes between each voxel and the voxels marked by the samples of the class to produce a minimum distance map that is localized. This map will have low values near the voxels of the area of interest (near in terms of both distance and intensity difference) and high values where it is computed to be outside the area of interest.

<a name="what-can-it-do"></a>
## What can it do?

Produce masks for the different areas specified on the sample of labels.

For brain tumor specifically, produce masks for:
- necrotic and non-enhancing tumor core
- GD-enhancing tumor 
- peritumoral edema
- other areas

See [advanced usage](EXTRA_ADVANCED.md) section for complementary operations.

See [alternative modes](EXTRA_MODES.md) section for alternative methods for segmentation.

<a name="how-to-use"></a>
## How to use?

### Input files

Supports 2D/3D NIfTI and DICOM images (.nii, .nii.gz, .dcm, .dicom files).

##### Input Images

The images can be MR, CT, PET or anything else (or a combination of all). The restriction is that all images should have exactly the same size and be co-registered images of one subject taken during a single session. 

<a name="#input-labels"></a>
##### Input Labels

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

See [here](EXTRA_ADVANCED.md#use-different-labels) for instructions on how to use different labels for the tumor-specific classes. 

Consider that the quality of the segmentations are highly correlated to the quality of the provided "sample of labels" image. Also, especially in the more difficult cases, there might be misclassified regions in the output segmentation. If that is the case the user can add more labels in the "sample of labels" images that are within these regions and run the tool again. 

[CaPTk](https://www.med.upenn.edu/cbica/captk/) can be used to draw the labels image. Drag and drop the MRI images in the "Images" area, draw labels through the "Drawing" tab and export using File>Save>ROI. For drawing more labels to an old labels image use File>Load>ROI, draw the new labels and export.

<a name="#run"></a>
### Run

###### Brain tumor

For brain tumor cases you can still use the general run instructions below, but better results can be achieved if:
- There is at least one of the following MRI modalities: FLAIR, T1, T1CE (Gd) or T2
- You supply the special modalities through their respective parameter
- You use the brain tumor specific labels [specified above](#brain-tumor)

In that case the usage, from within the same directory as the executable, is:

```console
GeodesicTraining.exe -flair "FULL_PATH_TO_FLAIR.nii.gz" -t1 "FULL_PATH_TO_T1.nii.gz" -t1ce "FULL_PATH_TO_T1CE.nii.gz" -t2 "FULL_PATH_TO_T2.nii.gz" -i "FULL_PATH_TO_EXTRA_IMAGE1.nii.gz,FULL_PATH_TO_EXTRA_IMAGE2.nii.gz,..." -l "FULL_PATH_TO_SAMPLE_OF_LABELS_IMAGE.nii.gz" -o "FULL_PATH_TO_DESIRED_OUTPUT_FOLDER"
```

- If any of FLAIR/T1/T1CE/T2 images are not available, don't use their respective parameter
- If no extra images are available, don't use the ```-i``` parameter
- "labels_res.nii.gz" will be the output segmentation with healthy tissue labels
- "labels_res_renamed.nii.gz" will be the output segmentation with healthy tissue labels removed
- In a UNIX environment replace "GeodesicTraining.exe" with "./GeodesicTraining".
- If the images are 2D use ```-id 2``` (although it might work without it).

###### General

From within the same directory as the executable:

```console
GeodesicTraining.exe -i "FULL_PATH_TO_IMAGE1.nii.gz,FULL_PATH_TO_IMAGE2.nii.gz,..." -l "FULL_PATH_TO_SAMPLE_OF_LABELS_IMAGE.nii.gz" -o "FULL_PATH_TO_DESIRED_OUTPUT_FOLDER"
```

- "labels_res.nii.gz" will be the output segmentation
- "labels_res_renamed.nii.gz" will be the output segmentation with healthy tissue labels removed __if__ you supply ```-cl L/0``` where L is the non-zero integer label used for healthy tissue.
- In a UNIX environment replace "GeodesicTraining.exe" with "./GeodesicTraining".
- If the images are 2D use ```-id 2``` (although it might work without it).

###### Advanced

For more advanced usage (comparing segmentation to ground truth, optional parameters, using different labels, isolating areas in an image and custom SVM configurations) see the [advanced usage](EXTRA_ADVANCED.md) section.

<a name="how-to-build"></a>
## How to build

Depends on ITK, OpenCV. Requires a C++11 compliant compiler and cmake. You have to pass the paths to ITK and OpenCV include directories to cmake. 

###### Linux

An example cmake command on linux could look like this (but the CMAKE_INSTALL_PREFIX and the ITK version might be different):

```console
cmake -DCMAKE_INSTALL_PREFIX="/usr/local" -DCMAKE_PREFIX_PATH="${CMAKE_INSTALL_PREFIX}/include/ITK-4.13;${CMAKE_INSTALL_PREFIX}/include/opencv4" -DCMAKE_MODULE_PATH="${CMAKE_INSTALL_PREFIX}/include/opencv4" ..
```

For more information on installing the dependencies on linux you can see [here for ITK](https://itk.org/Wiki/ITK_Configuring_and_Building_for_Ubuntu_Linux) and [here for OpenCV](https://docs.opencv.org/3.4/d7/d9f/tutorial_linux_install.html), but you can also try to look for newer versions.

###### Windows

This library can be built on Windows by providing the paths to ITK and OpenCV. See [here for ITK](https://itk.org/ITK/resources/software.html) and [here for OpenCV](https://docs.opencv.org/2.4/doc/tutorials/introduction/windows_install/windows_install.html), but you can also try to look for newer versions.

In general, the ```CMAKE_PREFIX_PATH``` should look something like this ```"C:/Libraries/InsightToolkit-4.13.1/build;C:/Libraries/opencv/build/x64/vc15/lib"``` and ```CMAKE_MODULE_PATH``` something like this ```"C:/Libraries/InsightToolkit-4.13.1/build/lib/cmake/ITK-4.13"```

If when trying to run the application you get prompted that ```opencv_world343.dll``` or  ```opencv_world343d.dll``` is missing, manually copy it from opencv's build path to the directory where the 'GeodesicTraining' executable is located. (The path to the dll will look something like this ```C:\Libraries\opencv\build\x64\vc15\bin```)

<a name="#how-to-use-as-a-library"></a>
## How to use as a library

The directory GeodesicTrainingSegmention of this repository can be used as a library. [See here for details](GeodesicTrainingSegmentation/README.md).
