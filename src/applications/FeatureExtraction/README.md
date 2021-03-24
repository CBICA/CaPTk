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

## Preprocessing

CaPTk's Feature Extraction provides default resampling to ensure all computations happen in the same physical space.

The ```Generic``` options in the parameter file control all options, which are:

| Parameter | Description |
|------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Resampling | Resamples all images and masks to this value of voxel/pixel resolution before computations (0: no resampling): reduce this number of more accuracy |
| ResamplingInterpolator_Image | Type of interpolator to use for input image(s) if resampling is happening; ignored if m_resamplingResolution = 0 |
| ResamplingInterpolator_Mask | Type of interpolator to use for the mask image if resampling is happening; ignored if m_resamplingResolution = 0 |
| Quantization_Extent | Whether the quantization of Intensities is supposed to happen on a per-image basis or on a per-ROI basis |
| Quantization_Type | FixedBinNumber (FBN): the bins are uniformly distributed across the minimum and maximum extent in the ROI/Image as defined under 'Quantization_Extent';<br> FixedBinSize (FBS): bins are added in a fixed step between the minimum and maximum extent in the ROI/Image as defined under 'Quantization_Extent' the requested size is provided in 'Bins';<br> Equal: each bin holds an equal number of intensities |
| SliceComputation | Controls whether non-Intensity features are calculated along the slice with the largest area along the 3 axes: valid for 3D images only |
| NaN-Handling | Specify how to handle features with NaN values 'Remove' removes the feature from output - keep in mind this might cause issues with output file in multi-subject (i.e. Training/Batch) mode |

### Making your own parameter file

1. Copy the sample parameter file
2. Remove the rows of features that you don't need (except intensity features, which are always calculated)
3. ```Column C```: Type, denotes what type of value it expects for the parameter
4. ```Column D```: Range, lists of possible values for the parameter
5. ```Column E```: Value of the specified parameter

### Effect of Resampling

Please see the following presentation: https://upenn.box.com/v/spacingsIssue

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


## Feature Extraction on the Image Processing Portal

You may also use Feature Extraction from the CBICA Image Processing Portal ("IPP"), accessible [here](https://ipp.cbica.upenn.edu/) .

The Multi-subject Feature Extraction experiment allows you to upload batches of images and ROIs. The experiment will automatically detect subject groupings, based on the filenames of the uploaded images and masks.

For example, if you upload AAAA_01012021_T2_LPS.nii.gz and AAAA_01012021_T1CE_LPS.nii.gz as images, and AAAA_01012021_mask_LPS.nii.gz as an ROI, those images will automatically be grouped in the output, using the same ROI for each modality. You may include multiple subjects with separate ROIs; subjects with different names will not interfere with one another. Even BBBB_T1_01012020.nii.gz and BBBB_T1_02022020.nii.gz will be detected as different subjects and will require separate ROI masks.
You must also specify ROI label values and names as you would with the Feature Extraction command-line interface.

This approach uses one single parameter file that is applied to all uploaded images/subjects. Please see the [sample parameter file](https://github.com/CBICA/CaPTk/blob/master/src/applications/FeatureExtraction/data/1_params_default.csv) .

Modality detection in the filenames relies on a known list of modality names. We detect some modalities by default (see the IPP page for details), but you can also specify custom modality names that will be detected in the "Filename Parsing" section.

## Frequently Asked Questions

- What do the various parameters mean and how do I customize them?

Take a look at the [sample parameter file](https://github.com/CBICA/CaPTk/blob/master/src/applications/FeatureExtraction/data/1_params_default.csv) for details.

- Why is preprocessing required?

See [this presentation](https://upenn.box.com/v/spacingsIssue) for a detailed explanation.

- Why do I see `#NAME` for some features?

These would be NaNs and depending on how features have been extracted, these are completely normal.

- Anything else?

Please [open a new issue](https://github.com/CBICA/CaPTk/issues/new?assignees=&labels=&template=bug-report.md&title=) 