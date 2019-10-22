cwlVersion: v1.0
class: CommandLineTool
baseCommand: WhiteStripe
inputs:
  inputImages:
    type: File
    label: 3D NIfTI Image(s)
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Images on which WhiteStripe needs to be applied.Delieanted by ','.
  output:
    type: string
    label: Dir with write access
    inputBinding:
      position: 1
      prefix: -o
    doc: Output directory where results are to be saved.Can be absolute path to file if input is a single image.
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
  radius:
    type: float?
    label: 0.0 to 5.0
    inputBinding:
      position: 1
      prefix: -r
    doc: "WhiteStripe Radius.Default: 0.05."
  skullStrippedImg:
    type: boolean?
    label: 1 or 0
    inputBinding:
      position: 1
      prefix: -sk
    doc: "Whether skull stripped image is passed.If not, 'zSliceRange' is needed.Default: 1."
  zSliceRange:
    type: int?
    label: 50 to 150
    inputBinding:
      position: 1
      prefix: -z
    doc: "z-slice Range for cropping.Delieanted by '-'.Default (if image is not skull-stripped): 80-120."
  tissuesMax:
    type: int?
    label: 5|10|20
    inputBinding:
      position: 1
      prefix: -t
    doc: "Max Tissues (5 or 10 or 20).Default: 5."
  maxSmooth:
    type: float?
    label: 0.0 to 50.0
    inputBinding:
      position: 1
      prefix: -m
    doc: "Max Smoothing.Default: 10.0."
  deltaSmooth:
    type: float?
    label: 0.0 to 10.0
    inputBinding:
      position: 1
      prefix: -d
    doc: "Smoothing Delta.Default: 0.5."
  binsHist:
    type: int?
    label: 100 to 3000
    inputBinding:
      position: 1
      prefix: -b
    doc: "Number of Histogram bins to do processing.Default: 2000."
  t1Image:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -t1
    doc: "T1 Image being passed or not.Default: 1."
hints:
  SoftwareRequirement:
    packages:
      WhiteStripe:
        version:
          - 1.7.3.nonRelease.20190819