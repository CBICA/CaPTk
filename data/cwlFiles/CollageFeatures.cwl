cwlVersion: v1.0
class: CommandLineTool
baseCommand: CollageFeatures
inputs:
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: Path to an input image from which features will be extracted..
  mask:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: Path to a mask that will be considered as binary. The highest pixel value will be considered as information and all other values will be considered outside the mask..
  outputfile:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: Path to the output CSV file..
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
  dimensions:
    type: int?
    label: 2-3
    inputBinding:
      position: 1
      prefix: -d
    doc: Optional number of dimensions upon which to run collage. Supported values are 2 and 3. If left out, we will default to the dimensionality of the image itself, which may not reflect expected behavior if the image has an alpha channel..
  svdradius:
    type: int?
    label: none
    inputBinding:
      position: 1
      prefix: -s
    doc: SVD radius is used for the dominant angle calculation pixel radius. DEFAULTS to 5 and is suggested to remain at the default..
  binsize:
    type: int?
    label: none
    inputBinding:
      position: 1
      prefix: -b
    doc: Number of bins to use while calculating the grey level cooccurence matrix. DEFAULTS to 64..
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
      CollageFeatures:
        version:
          - 1.9.0.Alpha