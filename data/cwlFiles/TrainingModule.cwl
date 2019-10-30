cwlVersion: v1.0
class: CommandLineTool
baseCommand: TrainingModule
inputs:
  features:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -f
    doc: The input file having features (*.csv)..
  label:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -l
    doc: The input file having target labels (*.csv)..
  classifier:
    type: int
    label: none
    inputBinding:
      position: 1
      prefix: -c
    doc: The SVM kernel to be used in developing model (1=Linear, 2=RBF)..
  configuration:
    type: int
    label: none
    inputBinding:
      position: 1
      prefix: -n
    doc: The Configuration type, Cross-validation (n=1), Split Train-Test (n=2), Train only (n=3), and Test only (n=4)..
  configuration parameters:
    type: int
    label: none
    inputBinding:
      position: 1
      prefix: -k
    doc: The number of folds for Cross-validation (5/10) and the size of training set for TrainTest (k<n)..
  output:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: The model direcory (needed only when n=4).
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
      TrainingModule:
        version:
          - 1.7.3.nonRelease.20190819