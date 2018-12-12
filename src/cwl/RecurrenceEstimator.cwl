cwlVersion: v1.0
class: CommandLineTool
version: 1.6.1
baseCommand: RecurrenceEstimator
inputs:
  type:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -t
    doc: The option of preparing a new model (=0), and for testing on an existing model (=1).
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: The input directory having test subjects.
  output:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: The output direcory to write output.
  model:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: The directory having SVM models.
  Logger:
    type: string?
    label: log file which user has write access to
    inputBinding:
      position: 1
      prefix: -L
    doc: Full path to log file to store console outputs.By default, only console output is generated.