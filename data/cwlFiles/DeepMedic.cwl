cwlVersion: v1.0
class: CommandLineTool
baseCommand: DeepMedic
inputs:
  T1CE:
    type: File
    label: none
    inputBinding:
      position: 1
      prefix: -t1c
    doc: The input T1CE or T1Gd image file..
  T1:
    type: File
    label: none
    inputBinding:
      position: 1
      prefix: -t1
    doc: The input T1 image file..
  FLAIR:
    type: File
    label: none
    inputBinding:
      position: 1
      prefix: -t2
    doc: The input T2 image file..
  output:
    type: Directory
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: The output File..
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
  mask:
    type: File?
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: The Optional input mask file..This is needed for normalization only.
  modelDir:
    type: Directory?
    label: none
    inputBinding:
      position: 1
      prefix: -md
    doc: The trained model to use.Defaults to 'CaPTk_installDir/data/deepMedic/brainSegmentation'.
  quantLower:
    type: float?
    label: 0-100
    inputBinding:
      position: 1
      prefix: -ql
    doc: "The Lower Quantile range to remove.This is needed for normalization only.Default: 5."
  quantUpper:
    type: float?
    label: 0-100
    inputBinding:
      position: 1
      prefix: -qu
    doc: "The Upper Quantile range to remove.This is needed for normalization only.Default: 95."
  cutOffLower:
    type: float?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -cl
    doc: "The Lower Cut-off (multiple of stdDev) to remove.This is needed for normalization only.Default: 3."
  cutOffUpper:
    type: float?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -cu
    doc: "The Upper Cut-off (multiple of stdDev) to remove.This is needed for normalization only.Default: 3."
  Logger:
    type: string?
    label: log file which user has write access to
    inputBinding:
      position: 1
      prefix: -L
    doc: Full path to log file to store console outputs.By default, only console output is generated.
hints:
  SoftwareRequirement:
    packages:
      DeepMedic:
        version:
          - 1.7.3.nonRelease.20190819