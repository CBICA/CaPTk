cwlVersion: v1.0
class: CommandLineTool
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
      PopulationAtlases:
        version:
          - 1.7.3.nonRelease.20190819