cwlVersion: v1.0
class: CommandLineTool
baseCommand: FeatureExtraction
inputs:
  outputDir:
    type: Directory
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: Absolute path of directory to save results.Result can be a CSV or Feature Maps (for lattice).
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
  paramFile:
    type: File?
    label: .csv
    inputBinding:
      position: 1
      prefix: -p
    doc: "A csv file with all features and its parameters filled.Default: '../data/1_params_default.csv'."
  batchFile:
    type: File?
    label: .csv
    inputBinding:
      position: 1
      prefix: -b
    doc: Input file with Multi-Patient Multi-Modality details.Header format is as follows:.'PATIENT_ID,IMAGES,MASK,ROI,SELECTED_ROI,ROI_LABEL,.SELECTED_FEATURES,PARAM_FILE'.Delineate individual fields by '|'.
  name_patient:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -n
    doc: Patient id.Required for single subject mode.
  imagePaths:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: "Absolute path of each coregistered modality.Delineate by ','.Example: -i c:/test1.nii.gz,c:/test2.nii.gz.Required for single subject mode."
  imageTypes:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -t
    doc: "Names of modalities to be processed.Delineate by ','.Example: -t T1,T1Gd.Required for single subject mode."
  maskPath:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: Absolute path of mask coregistered with images.Required for single subject mode.
  roi:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -r
    doc: "List of roi for which feature extraction is to be done.Delineate by ','.Example: -r 1,2.Required for single subject mode."
  labels:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -l
    doc: "Labels variables for selected roi numbers.Delineate by ','.Usage: -l Edema,Necrosis.Required for single subject mode."
  verticalConc:
    type: boolean?
    label: flag
    inputBinding:
      position: 1
      prefix: -vc
    doc: Whether vertical concatenation is needed or not.Horizontal concatenation is useful for training.Defaults to '0'.
  featureMapsWrite:
    type: boolean?
    label: flag
    inputBinding:
      position: 1
      prefix: -f
    doc: Whether downsampled feature maps are written or not.For Lattice computation ONLY.Defaults to '0'.
  threads:
    type: int?
    label: 1-64
    inputBinding:
      position: 1
      prefix: -th
    doc: Number of (OpenMP) threads to run FE on.Defaults to '1'.This gets disabled when lattice is disabled.
  offsets:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -of
    doc: "Exact offset values to pass on for GLCM & GLRLM.Should be same as ImageDimension and in the format '<offset1>,<offset2>,<offset3>'.This is scaled on the basis of the radius.Example: '-of 0x0x1,0x1x0'."
  debug:
    type: boolean?
    label: True or False
    inputBinding:
      position: 1
      prefix: -d
    doc: Whether to print out additional debugging info.Defaults to '0'.
  debugWrite:
    type: boolean?
    label: True or False
    inputBinding:
      position: 1
      prefix: -dw
    doc: Whether to write intermediate files or not.Defaults to '0'.
  Logger:
    type: File?
    label: Text file with write access
    inputBinding:
      position: 1
      prefix: -L
    doc: Full path to log file to store logging information.By default, only console output is generated.
  unitTest:
    type: File?
    label: Path to reference output
    inputBinding:
      position: 1
      prefix: -ut
    doc: Whether to run unit test or not.Disabled for batch processing.
hints:
  SoftwareRequirement:
    packages:
      FeatureExtraction:
        version:
          - 1.7.3.nonRelease.20190819