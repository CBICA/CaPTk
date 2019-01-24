# Feature Extraction

## Command Line Usage

### The Interface


### Windows PowerShell / CMD:
```
FeatureExtraction.exe –o C:\path\to\output\csv –b C:\path\to\batch\csv –p C:\path\to\parameters\csv [– d 1] [– vc 1]
```

### Linux/macOS Terminal:
```
FeatureExtraction –o /path/to/output/csv –b /path/to/batch/csv –p /path/to/parameters/csv [– d 1] [– vc 1]
```

## Locating the executable

### Windows

Under the CaPTk installation tree, you will see the folder */bin/*; the FeatureExtraction executable will be there.

### macOS

The macOS package is an [Bundle](https://en.wikipedia.org/wiki/Bundle_(macOS)#macOS_application_bundles) and the FeatureExtraction executable can be access using the following [cwl](https://www.commonwl.org/)-enabled command (assuming you have the CaPTk AppImage in your PATH): ```CaPTk.app featureextraction $your_command_goes_here```

### Linux

The Linux package is an [AppImage](https://appimage.org/) and the FeatureExtraction executable can be access using the following [cwl](https://www.commonwl.org/)-enabled command (assuming you have the CaPTk AppImage in your PATH): ```captk featureextraction $your_command_goes_here```

## The Parameter File

Controls which features to calculate its respective parameters.

Sample: ```$CaPTk_Install_Directory/data/features/1_params_default.csv```

### Making your own parameter file

1. Copy the sample parameter file
2. Remove the rows of features that you don't need (Except intensity in row 2, which is always calculated)
3. Column C: Type, denotes what type of value it expects for the parameter
4. Column D: Range, lists of possible values for the parameter
5. Column E: Value of the specified parameter

## Batch Processing

Calculating features (using the same parameter file) for a list of subjects (with multi-modal data).

Sample: ```$CaPTk_Install_Directory/share/featureExtractionBatch/batch_extraction.csv```

1. Column A: Patient, is the naming to ID the subject
2. Column B: Images, is the full path(s) to images, separated by "|" as delimiter
3. Column C: Modalities, is the naming for each images, separated by "|" as delimiter
4. Column D: ROIFile, is the full path to a label file (binary of multi-label) for detecting ROI
5. Column E: SELECTED_ROI, is the value(s) in the label file that you want to use as ROI region, separated by "|" as delimiter
6. Column F: ROI, is the corresponding naming of the label(s) you chose in Column E
7. [OPTIONAL] Column G: OUTPUT, the optional output file, if you want each subject's output to be written into a different file