cwlVersion: v1.0
class: CommandLineTool
version: 1.6.2.Beta
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
    doc: The Configuration type, either Cross-validation or Split Train-Test (1=CrossValidtion, 2=TrainTest)..
  configuration parameters:
    type: int
    label: none
    inputBinding:
      position: 1
      prefix: -k
    doc: The number of folds for Crossvalidation (5/10) and the size of training set for TrainTest (k<n)..
  output:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: The output direcory to write output.
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