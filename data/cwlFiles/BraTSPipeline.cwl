cwlVersion: v1.0
class: CommandLineTool
baseCommand: BraTSPipeline
inputs:
  t1ceImage:
    type: string
    label: Input Image (DICOM or NIfTI)
    inputBinding:
      position: 1
      prefix: -t1c
    doc: Input structural T1-weighted post-contrast image.
  t1Image:
    type: string
    label: Input Image (DICOM or NIfTI)
    inputBinding:
      position: 1
      prefix: -t1
    doc: Input structural T1-weighted pre-contrast image.
  t2Image:
    type: string
    label: Input Image (DICOM or NIfTI)
    inputBinding:
      position: 1
      prefix: -t2
    doc: Input structural T2-weighted contrast image.
  flImage:
    type: string
    label: Input Image (DICOM or NIfTI)
    inputBinding:
      position: 1
      prefix: -fl
    doc: Input structural FLAIR contrast image.
  outputDir:
    type: Directory
    label: Directory
    inputBinding:
      position: 1
      prefix: -o
    doc: Output directory for final output.
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
  skullStrip:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -s
    doc: "Flag whether to skull strip or not.Defaults to 1.This uses DeepMedic: https://cbica.github.io/CaPTk/seg_DL.html."
  brainTumor:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -b
    doc: "Flag whether to segment brain tumors or not.Defaults to 1.This uses DeepMedic: https://cbica.github.io/CaPTk/seg_DL.html."
  debug:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -d
    doc: Print debugging information.Defaults to 1.
  interFiles:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -i
    doc: Save intermediate files.Defaults to 1.
hints:
  SoftwareRequirement:
    packages:
      BraTSPipeline:
        version:
          - 1.8.0.Beta