cwlVersion: v1.0
class: CommandLineTool
baseCommand: GreedyRegistration
inputs:
  movingImage:
    type: File
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Image for processing.Becomes moving image in Registration mode.
  fixedImage:
    type: string
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -f
    doc: Fixed Image for registration.
  outputImage:
    type: File
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -o
    doc: Output Image for processing.
  matrix:
    type: string
    label: N.A
    inputBinding:
      position: 1
      prefix: -t
    doc: Registration Matrix.
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
  registration:
    type: string?
    label: N.A
    inputBinding:
      position: 1
      prefix: -reg
    doc: Switch to registration mode.
  transformation:
    type: string?
    label: N.A
    inputBinding:
      position: 1
      prefix: -trf
    doc: Switch to transformation mode.
  affine:
    type: string?
    label: N.A
    inputBinding:
      position: 1
      prefix: -a
    doc: Affine Registration(Default).
  rigid:
    type: string?
    label: N.A
    inputBinding:
      position: 1
      prefix: -r
    doc: Rigid Registration.
  metrics:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: "MI: mutual information.NMI(Default): normalized mutual information.NCC -r 2x2x2: normalized cross-correlation."
  radius:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -ri
    doc: "Patch radius for metrics.Eg: 2x2x2."
  greedyIterations:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -n
    doc: "Number of iterations per level of multi-res (Default: 100x50x5).Corresponds to low level, Mid Level and High Level resolution.Pattern: NxNxN."
  threads:
    type: int?
    label: none
    inputBinding:
      position: 1
      prefix: -th
    doc: Number of threads for algorithm.If not suppllied gets set to default 4.
hints:
  SoftwareRequirement:
    packages:
      GreedyRegistration:
        version:
          - 1.7.3.nonRelease.20190819