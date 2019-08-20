cwlVersion: v1.0
class: CommandLineTool
baseCommand: EGFRvIIISurrogateIndex
inputs:
  image:
    type: File
    label: NIfTI or DICOM
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Perfusion image on which computation is done.
  mask:
    type: File
    label: NIfTI or DICOM
    inputBinding:
      position: 1
      prefix: -m
    doc: Mask containing near (1) and far (2) labels.
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
hints:
  SoftwareRequirement:
    packages:
      EGFRvIIISurrogateIndex:
        version:
          - 1.7.3.nonRelease.20190819