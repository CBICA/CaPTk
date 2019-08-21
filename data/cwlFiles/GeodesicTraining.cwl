cwlVersion: v1.0
class: CommandLineTool
baseCommand: GeodesicTraining
inputs:
  runtest:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -rt
    doc: Runs the tests.
  cwl:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -cwl
    doc: Generates a .cwl file for the software.
  mode:
    type: string?
    label: text
    inputBinding:
      position: 1
      prefix: -m
    doc: "Available modes are: \n\n                             reversegeotrain -> [DEFAULT MODE] use AGD and SVMs to produce labels\n                             svmlabels       -> use SVMs to produce labels\n           [2-class]         svmpseudo       -> use SVMs to produce pseudoprobability maps\n           [1-class]         agd             -> run AGD to produce agd maps\n           [2-class]         geotrain        -> use SVMs and AGD to produce agd maps and threshold them\n           [2-class]         geotrainfull    -> AGD ---> SVMs ---> AGD ---> threshold\n                             rf              -> Random Forest\n                             rfauto          -> Random Forest with automatic parameter tuning\n                             agdrf           -> Random Forest with AGD input images\n                             agdrfauto       -> Random Forest with AGD input images and automatic parameter tuning\n           [util]            segment         -> segment an area from an image according to a labels image\n           [util]            labelthres      -> convert AGD maps to labels images\n           [util]            checkaccuracy   -> Check the accuracy of a segmentation in relation to a ground truth one\n           [util]            changelabels    -> change the labels of a labels image\n           [util]            generateconfig  -> generate a configuration file from pretrained svm models\n."
  input:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -i
    doc: List of full paths to the input files, separated by comma..
  Logger:
    type: string?
    label: log file with user write access
    inputBinding:
      position: 1
      prefix: -L
    doc: Full path to log file to store console output..By default, only console output is generated..
  labels:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -l
    doc: Full path to the input labels file..
  output:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -o
    doc: Directory for output..
  configuration:
    type: string?
    label: .yaml
    inputBinding:
      position: 1
      prefix: -c
    doc: Full path to SVM Configuration file..
  rfconfiguration:
    type: string?
    label: .config
    inputBinding:
      position: 1
      prefix: -rfc
    doc: Full path to RF Configuration file..
  datasetname:
    type: string?
    label: text
    inputBinding:
      position: 1
      prefix: -d
    doc: The dataset name.Used for better output organization..Folders will be created inside the output dir..
  tag:
    type: string?
    label: text
    inputBinding:
      position: 1
      prefix: -tg
    doc: A tag for this execution.Used for better output organization..Tag will be used only if dataset name is set..
  nodatetime:
    type: string?
    label: "-"
    inputBinding:
      position: 1
      prefix: -nd
    doc: Don't include datetime in the directory name.Used for better output organization..Datetime will still be used if datasetname is provided but no tag..
  noclose:
    type: string?
    label: "-"
    inputBinding:
      position: 1
      prefix: -n
    doc: Option not to close the window after the program has finished..
  reportseconds:
    type: string?
    label: "-"
    inputBinding:
      position: 1
      prefix: -r
    doc: Enable time report..Creates time_report.txt file in output directory..
  threshold:
    type: string?
    label: float
    inputBinding:
      position: 1
      prefix: -t
    doc: Threshold value.
  changelabelsmap:
    type: string?
    label: list of tuples
    inputBinding:
      position: 1
      prefix: -cl
    doc: "Importance for each model separated by comma..Example: Change label 1 -> 2, and labels 3 and 4 -> 5.the changelabelsmap value should be \"1/2,3/5,4/5\"."
  inputmodels:
    type: string?
    label: .xml
    inputBinding:
      position: 1
      prefix: -im
    doc: "List of input models separated by comma..For mode: \"generateconfig\"."
  importancevalue:
    type: string?
    label: double vector
    inputBinding:
      position: 1
      prefix: -iv
    doc: "Importance for each model separated by comma..For mode: \"generateconfig\"."
  groundtruth:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -g
    doc: Ground truth image (for checking accuracy)..
  groundtruthskip:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -gs
    doc: List of labels to skip when checking accuracy .in relation to ground truth image. Default is 0..
  threads:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -j
    doc: Number of threads that will be used..Give "max" for infinite..Default is 16..
  flair:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -flair
    doc: Flair image (For MRI images only).Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images..
  t1:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -t1
    doc: T1 image (For MRI images only).Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images..
  t1ce:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -t1ce
    doc: T1ce image (For MRI images only).Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images..
  t2:
    type: string?
    label: .nii.gz
    inputBinding:
      position: 1
      prefix: -t2
    doc: T2 image (For MRI images only).Do not supply images in the i parameter for flair/t1/t1ce/t2 MRI images..
  labeltc:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -ltc
    doc: Label for tumor core (For MRI images only).
  labelet:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -let
    doc: Label for enhanced tumor (For MRI images only).
  labeled:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -led
    doc: Label for edema (For MRI images only).
  labelht:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -lht
    doc: Label for healthy tissue (For MRI images only).
  labelofinterest:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -loi
    doc: Label of interest..Default is 1..Any other non-zero label can be used for areas that are not of interest..
  nosubsample:
    type: string?
    label: "-"
    inputBinding:
      position: 1
      prefix: -ns
    doc: Option to not subsample if the number of samples is high (SVM).
  maxsamples:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -ms
    doc: Max samples to use (SVM Subsampling).
  nobalancesamples:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -nb
    doc: Don't balance the subsampling (SVM Subsampling).
  imagedimensions:
    type: string?
    label: int
    inputBinding:
      position: 1
      prefix: -id
    doc: Input image(s) dimensions [only 3D supported for now].
hints:
  SoftwareRequirement:
    packages:
      GeodesicTraining:
        version:
          - 1.7.3.nonRelease.20190819