cwlVersion: v1.0
class: CommandLineTool
version: 1.6.1
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
  No. of folds:
    type: int
    label: none
    inputBinding:
      position: 1
      prefix: -k
    doc: The number of folds to develop the model (5/10)..
  output:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: The output direcory to write output.
  Logger:
    type: string?
    label: log file which user has write access to
    inputBinding:
      position: 1
      prefix: -L
    doc: Full path to log file to store console outputs.By default, only console output is generated.