# Advanced Usage

### Contents 
- [Use different labels](#use-different-labels)
- [Change labels in a labels image](#change-labels-in-a-labels-image)
- [Unknown modality](#unknown-modality)
- [Compare segmentation to ground truth](#compare-segmentation-to-ground-truth)
- [More optional parameters](#more-optional-parameters)
- [Isolate an area in an image](#isolate-an-area-in-an-image)
- [Custom SVM ensemble](#custom-svm-ensemble)

<a name="#use-different-labels"></a>
## Use different labels

If you want to use different labels than the ones specified [here](README.md#input-labels) use parameters:
* ```-ltc LABEL``` For specifying a label for tumor core.
* ```-let LABEL``` For specifying a label for enhancing tumor.
* ```-led LABEL``` For specifying a label for edema.
* ```-lht LABEL``` For specifying a label for healthy tissue. 

where LABEL is an non-zero integer.

<a name="#change-labels-in-a-labels-image"></a>
## Change labels in a labels image

```console
GeodesicTraining.exe -m changelabels -l "FULL_PATH_TO_LABELS_IMAGE.nii.gz" -o "FULL_PATH_TO_DESIRED_OUTPUT_FOLDER" -cl CHANGE_LABELS_STRING
```

where CHANGE_LABELS_STRING should specify the desired changes. For example if label 1 need to be changed to 2 and labels 3 and 4 need to be changed to 5 the CHANGE_LABELS_STRING should be: "1/2,3/5,4/5"

```-cl CHANGE_LABELS_STRING``` can also be used in the command described [here](README.md/#run). Note that the segmentation is always saved as "labels_res.nii.gz" and the segmentation with some labels changed is saved as "labels_res_renamed.nii.gz". The default CHANGE_LABELS_STRING is "X/0" where X is the label for healthy tissue (default label for HT is 3), so be sure to include to a custom CHANGE_LABELS_STRING if removing healthy brain tissue labels is still desired to be performed.

<a name="#unknown-modality"></a>
## Unknown modality

The modalities with specific support are T1, T2, T1CE (Gd) and FLAIR. Images provided using these modalities should be supplied using the parameters ```-t1```, ```-t2```, ```-t1ce``` or ```-flair``` respectively, because different operations happen on different modalities. If there are more images available using other modalities that might be useful they can be provided as input using the parameter ```-i``` separated by comma like this:

```
-i "FULL_PATH_TO_EXTRA_IMAGE1,FULL_PATH_TO_EXTRA_IMAGE2,..."
```

Please note that using ```-i``` to supply T1, T2, T1CE or FLAIR will result in worst segmentations than when supplying them through the correct parameter.

<a name="#compare-segmentation-to-ground-truth"></a>
## Compare segmentation to ground truth

Please not that label 0 is the null label and is not used at the calculations.

```console
GeodesicTraining.exe -m checkaccuracy -l "FULL_PATH_TO_LABELS_IMAGE.nii.gz" -g "FULL_PATH_TO_GROUND_TRUTH_LABELS_IMAGE" -o "FULL_PATH_TO_DESIRED_OUTPUT_FOLDER"
```

Optional parameters for comparing to ground truth:
* ```-cl CHANGE_LABELS_STRING``` If you want to transform some of the labels of the input labels image to something else. For example changing label 1 to 2 and labels 3 and 4 to 5 is done using CHANGE_LABELS_STRING: "1/2,3/5,4/5". Please note that this parameter doesn't change anything to the ground truth image.
* ```-gs SKIP_LABELS_STRING``` If you want to *not* consider some labels when comparing. Label 0 is not considered regardless.

<a name="#more-optional-parameters"></a>
## More optional parameters

* ```--reportseconds``` Records duration and saves it to "time_report.txt"
* ```-j NUMBER_OF_THREADS``` Sets the number of threads manually. Default is 16.
* ```--noclose``` Option not to close the window after the program has finished.
* ```-d DATASET_NAME``` Useful for a better output folder organization. Using this parameter will create a folder inside output folder with the name specified with variable DATASET_NAME. Inside this directory another folder will be created every time the program is run with current datetime as name and it will contain the output of the execution. Alongside ```-d DATASET_NAME``` the following parameters can be used:
  * ```-tg TAG``` Alongside datetime what is specified in TAG will be part of the directory name.
  * ```--nodatetime``` Don't use datetime and use only tag for directory name. Will be ignored if no tag is provided.

<a name="#isolate-an-area-in-an-image"></a>
## Isolate an area in an image

```console
GeodesicTraining.exe -m segment -i "FULL_PATH_TO_INPUT_IMAGE" -l "FULL_PATH_TO_LABELS_IMAGE.nii.gz" -ll LABEL -o "FULL_PATH_TO_DESIRED_OUTPUT_FOLDER"
```

The output will be the input image where the voxels that *don't* have the label specified in variable LABEL become zero.

<a name="#custom-svm-ensemble"></a>
## Custom SVM ensemble

By default no configuration file need to be provided and an ensemble of SVMs with RBF, chi2 and histogram intersection kernels are trained, but a custom YAML configuration file can be provided using ```-c "FULL_PATH_TO_SVM_CONF_FILE"```

For instruction on what the .yaml file should contain, [see here](GeodesicTrainingSegmentation/EXTRA_ADVANCED.md/#the-different-keywords-are).
