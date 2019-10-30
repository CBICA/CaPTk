cwlVersion: v1.0
class: CommandLineTool
baseCommand: PerfusionPCA
inputs:
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: The input directory..
  type:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -t
    doc: The option of preparing a new model (=0), and for testing on an existing model (=1).
  number of PCAs:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -n
    doc: The number of principal components..
  output:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: The output directory..
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
      PerfusionPCA:
        version:
          - 1.7.3.nonRelease.20190819