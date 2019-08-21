cwlVersion: v1.0
class: CommandLineTool
baseCommand: SBRT_LungField
inputs:
  petImage:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -p
    doc: Absolute path of PET image.For example /data/.../pet.nii.gz.
  ctImage:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -c
    doc: Absolute path of CT image.For example /data/.../ct.nii.gz.
  outputImage:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: Absolute path and basename of output file (without extension).For example /output/.../label.
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
  maskImage:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: Absolute path of mask image.For example /data/.../mask.nii.gz.
  logFile:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -L
    doc: Absolute path of log file.For example log_file.txt.
hints:
  SoftwareRequirement:
    packages:
      SBRT_LungField:
        version:
          - 1.7.3.nonRelease.20190819