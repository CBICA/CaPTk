cwlVersion: v1.0
class: CommandLineTool
baseCommand: SBRT_Analysis
inputs:
  inputImage:
    type: File
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: Absolute path of PET image.For example pet.nii.gz.
  maskImage:
    type: File
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: Absolute path of mask image.For example mask.nii.gz.
  label:
    type: int
    label: none
    inputBinding:
      position: 1
      prefix: -l
    doc: Label value of the ROI.For example 2.
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
  outputFile:
    type: File?
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: Absolute path and basename of output file (without extension).For example radiomic_feature.
  logFile:
    type: File?
    label: none
    inputBinding:
      position: 1
      prefix: -L
    doc: Absolute path of log file.For example log_file.txt.
  Directory:
    type: Directory?
    label: none
    inputBinding:
      position: 1
      prefix: -D
    doc: Absolute path of model directory.For example C:/Model.
hints:
  SoftwareRequirement:
    packages:
      SBRT_Analysis:
        version:
          - 1.7.3.nonRelease.20190819