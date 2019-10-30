cwlVersion: v1.0
class: CommandLineTool
baseCommand: DirectionalityEstimate
inputs:
  labelMap:
    type: File
    label: NIfTI file
    inputBinding:
      position: 1
      prefix: -l
    doc: The label map on which calculations are to be done.
  outputFile:
    type: File
    label: text file
    inputBinding:
      position: 1
      prefix: -o
    doc: Where output gets stored.
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
  real:
    type: float?
    label: labelMap Dimensions
    inputBinding:
      position: 1
      prefix: -r
    doc: The real world coordinates (in mm) of file.Take example from tissue point file.Needs to be in same dimensionality as labelMap.Delineation is done using ','.
  index:
    type: string?
    label: labelMap Dimensions
    inputBinding:
      position: 1
      prefix: -i
    doc: The image indeces (in voxels) of file.Needs to be in same dimensionality as labelMap.Delineation is done using ','.
hints:
  SoftwareRequirement:
    packages:
      DirectionalityEstimate:
        version:
          - 1.7.3.nonRelease.20190819