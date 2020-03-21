# Feature Extraction

## Command Line Usage

The examples given are with the assumption that the user is running the application on Windows. Remove the *.exe* at the end for Linux/macOS.

### Example for Single Subject Computation

```
FeatureExtraction.exe \
-p C:/test/captk/params_test.csv \ # the parameter file which defines the features to extract and the related parameters
-n SubjectID \ # the subject ID of the current patient - used in the output csv
-i C:/test/captk/image.nii.gz \ # the input image(s) - multiple can be passed with ',' used as delimiter but need to be co-registered
-t FL \ # the modality of the input images (used in the output csv) - multiple can be passed with ',' used as delimiter
-m C:/test/captk/roi.nii.gz \ # the annotated mask coregistered with the input image(s)
-r 1 \ # the values in the annotated mask on which feature extraction is to occure -  multiple can be passed with ',' used as delimiter
-l TT \ # the labels of the selected ROIs (used in the output csv) - multiple can be passed with ',' used as delimiter
-vc 1 \ # vertically concatenate the output; if this is not passed, the output is horizontally concatenated to make it easier to read as a feature vector
-o C:/test/captk/output.csv # the output file
```

### Primary Help
```
FeatureExtraction.exe –u 
```

### Detailed Help
```
FeatureExtraction.exe –h
```

## Locating the executable

### Windows

Under the CaPTk installation tree, you will see the folder */bin/*; the FeatureExtraction executable will be there.

### macOS

The macOS package is a [Bundle](https://en.wikipedia.org/wiki/Bundle_(macOS)#macOS_application_bundles) and the FeatureExtraction executable can be access using the following [cwl](https://www.commonwl.org/)-enabled command (assuming you have the CaPTk AppImage in your PATH): ```CaPTk.app featureextraction $your_command_goes_here```

### Linux

The Linux package is an [AppImage](https://appimage.org/) and the FeatureExtraction executable can be access using the following [cwl](https://www.commonwl.org/)-enabled command (assuming you have the CaPTk AppImage in your PATH): ```captk featureextraction $your_command_goes_here```

## The Parameter File

Controls which features to calculate its respective parameters.

[Sample](https://github.com/CBICA/CaPTk/blob/master/src/applications/FeatureExtraction/data/1_params_default.csv): ```$CaPTk_Install_Directory/data/features/1_params_default.csv```

### Making your own parameter file

1. Copy the sample parameter file
2. Remove the rows of features that you don't need (except intensity features, which are always calculated)
3. ```Column C```: Type, denotes what type of value it expects for the parameter
4. ```Column D```: Range, lists of possible values for the parameter
5. ```Column E```: Value of the specified parameter

## Batch Processing

Calculating features (using the same parameter file) for a list of subjects (with multi-modal data).

[Sample](https://github.com/CBICA/CaPTk/blob/master/src/applications/FeatureExtraction/data/batchMode/batch_featureExtraction.csv): ```$CaPTk_Install_Directory/share/featureExtractionBatch/batch_extraction.csv```

1. ```Column A```: Patient, is the naming to ID the subject
2. ```Column B```: Images, is the full path(s) to images, separated by "|" as delimiter
3. ```Column C```: Modalities, is the naming for each images, separated by "|" as delimiter
4. ```Column D```: ROIFile, is the full path to a label file (binary of multi-label) for detecting ROI
5. ```Column E```: SELECTED_ROI, is the value(s) in the label file that you want to use as ROI region, separated by "|" as delimiter
6. ```Column F```: ROI, is the corresponding naming of the label(s) you chose in Column E
7. [OPTIONAL] ```Column G```: OUTPUT, the optional output file, if you want each subject's output to be written into a different file