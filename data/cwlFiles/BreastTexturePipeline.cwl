cwlVersion: v1.0
class: CommandLineTool
baseCommand: BreastTexturePipeline
inputs:
  inputImage:
    type: File
    label: DICOM
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Image for processing.
  outputDir:
    type: Directory
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -o
    doc: Dir with write access.All output files are written here.
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
  debugMode:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -d
    doc: "Enabled debug mode.Default: 0."
  resize:
    type: int?
    label: 0 - 100
    inputBinding:
      position: 1
      prefix: -r
    doc: "What resizing factor is to be applied.Default: 100."
hints:
  SoftwareRequirement:
    packages:
      BreastTexturePipeline:
        version:
          - 1.7.3.nonRelease.20190819