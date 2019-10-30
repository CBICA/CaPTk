cwlVersion: v1.0
class: CommandLineTool
baseCommand: GeodesicSegmentation
inputs:
  output:
    type: string
    label: NIfTI or Directory Path
    inputBinding:
      position: 1
      prefix: -o
    doc: Output File to save the results for single subject.If 'c' is passed, then this expects a directory where results will be stored.
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
  inputCSV:
    type: File?
    label: csv File
    inputBinding:
      position: 1
      prefix: -c
    doc: CSV containing list of input subjects and GLISTR/GLISTRboost labels.If this is provided, only 'o' is needed as parameter.
  imageCSV:
    type: string?
    label: String
    inputBinding:
      position: 1
      prefix: -cm
    doc: CSV header of mask correspoding to input image.Required if 'c' is passed.
  image:
    type: File?
    label: NIfTI or DICOM
    inputBinding:
      position: 1
      prefix: -i
    doc: Input image which needs to be segmented.
  mask:
    type: File?
    label: NIfTI or DICOM
    inputBinding:
      position: 1
      prefix: -m
    doc: Seed points in the image from where to start computation.
  normalize:
    type: boolean?
    label: flag
    inputBinding:
      position: 1
      prefix: -n
    doc: Normalize output map between 0-255 or not.Defaults to true/1.
  threshold:
    type: int?
    label: 0-255
    inputBinding:
      position: 1
      prefix: -t
    doc: Threshold distance for geodesic mask.By default, full geodesic mask is written.
hints:
  SoftwareRequirement:
    packages:
      GeodesicSegmentation:
        version:
          - 1.7.3.nonRelease.20190819