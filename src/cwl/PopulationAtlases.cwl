cwlVersion: v1.0
class: CommandLineTool
version: 1.6.1
baseCommand: PopulationAtlases
inputs:
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: The input directory..
  label:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -l
    doc: The input label file in .csv format..
  atlas:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -a
    doc: The atlas template..
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