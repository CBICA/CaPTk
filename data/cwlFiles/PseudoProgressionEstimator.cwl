cwlVersion: v1.0
class: CommandLineTool
baseCommand: PseudoProgressionEstimator
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
hints:
  SoftwareRequirement:
    packages:
      PseudoProgressionEstimator:
        version:
          - 1.7.3.nonRelease.20190819