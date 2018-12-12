cwlVersion: v1.0
class: CommandLineTool
version: 1.6.1
baseCommand: PerfusionPCA
inputs:
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: The input DSC-MRI image..
  mask:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: The input mask..
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
  Logger:
    type: string?
    label: log file which user has write access to
    inputBinding:
      position: 1
      prefix: -L
    doc: Full path to log file to store console outputs.By default, only console output is generated.