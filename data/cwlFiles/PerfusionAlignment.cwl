cwlVersion: v1.0
class: CommandLineTool
baseCommand: PerfusionAlignment
inputs:
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: The input DSC-MRI image..
  dicom file:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -a
    doc: The number of time-points after drop..
  t1ce file:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -c
    doc: The input T1 post-weighted image..
  echo time:
    type: float
    label: none
    inputBinding:
      position: 1
      prefix: -e
    doc: Echo time..
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
      PerfusionAlignment:
        version:
          - 1.7.3.nonRelease.20190819