cwlVersion: v1.0
class: CommandLineTool
baseCommand: DiffusionDerivatives
inputs:
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: The input DWI file..
  mask:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: The input mask file..
  Bval:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -b
    doc: The input bval file..
  Axial Diffusivity:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -a
    doc: The Axial Diffusivity image (1=YES, 0=NO, 1 (Default)).
  Fractional Anisotropy:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -f
    doc: The Fractional Anisotropy image (1=YES, 0=NO, 1 (Default)).
  Radial Diffusivity:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -r
    doc: The Radial Diffusivity image (1=YES, 0=NO, 1 (Default)).
  Apparent Diffusion Coefficient:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -t
    doc: The Apparent Diffusion Coefficient (1=YES, 0=NO, 1 (Default)).
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
      DiffusionDerivatives:
        version:
          - 1.7.3.nonRelease.20190819