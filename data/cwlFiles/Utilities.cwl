cwlVersion: v1.0
class: CommandLineTool
baseCommand: Utilities
inputs:
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
  inputImage:
    type: string?
    label: File or Dir
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Image (all CaPTk supported images) for processing.Directory to a single series DICOM only.
  maskImage:
    type: string?
    label: File or Dir
    inputBinding:
      position: 1
      prefix: -m
    doc: Input Mask (all CaPTk supported images) for processing.Directory to a single series DICOM only.
  outputImage:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -o
    doc: Output Image for processing.
  dicomDirectory:
    type: Directory?
    label: none
    inputBinding:
      position: 1
      prefix: -df
    doc: Absolute path of directory containing single dicom series.
  resize:
    type: int?
    label: 10-500
    inputBinding:
      position: 1
      prefix: -r
    doc: "Resize an image based on the resizing factor given.Example: -r 150 resizes inputImage by 150%.Defaults to 100, i.e., no resizing.Resampling can be done on image with 100."
  resizeResolution:
    type: float?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -rr
    doc: "[Resample] Isotropic resolution of the voxels/pixels to change to.Resize value needs to be 100.Defaults to 1.0."
  resizeInterp:
    type: string?
    label: NEAREST:LINEAR:BSPLINE:BICUBIC
    inputBinding:
      position: 1
      prefix: -ri
    doc: The interpolation type to use for resampling or resizing.Defaults to LINEAR.
  sanityCheck:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -s
    doc: Do sanity check of inputImage with the file provided in with this parameter.Performs checks on size, origin & spacing.Pass the target image after '-s'.
  information:
    type: boolean?
    label: true or false
    inputBinding:
      position: 1
      prefix: -inf
    doc: Output the information in inputImage.If DICOM file is detected, the tags are written out.
  cast:
    type: string?
    label: (u)char, (u)int, (u)long, (u)longlong, float, double
    inputBinding:
      position: 1
      prefix: -c
    doc: "Change the input image type.Examples: '-c uchar', '-c float', '-c longlong'."
  uniqueVals:
    type: boolean?
    label: true or false
    inputBinding:
      position: 1
      prefix: -un
    doc: Output the unique values in the inputImage.Pass value '1' for ascending sort or '0' for no sort.Defaults to '1'.
  boundingBox:
    type: File?
    label: NIfTI Mask
    inputBinding:
      position: 1
      prefix: -b
    doc: Extracts the smallest bounding box around the mask file.With respect to inputImage.Writes to outputImage.
  boundingIso:
    type: boolean?
    label: Isotropic Box or not
    inputBinding:
      position: 1
      prefix: -bi
    doc: Whether the bounding box is Isotropic or not.Defaults to true.
  testBase:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -tb
    doc: Baseline image to compare inputImage with.
  testRadius:
    type: int?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -tr
    doc: Maximum distance away to look for a matching pixel.Defaults to 0.
  testThresh:
    type: float?
    label: 0-5
    inputBinding:
      position: 1
      prefix: -tt
    doc: Minimum threshold for pixels to be different.Defaults to 0.0.
  createMask:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -cm
    doc: "Create a binary mask out of a provided (float) thresholds.Format: -cm lower,upper.Output is 1 if value >= lower or <= upper.Defaults to 1,Max."
  changeValue:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -cv
    doc: "Change the specified pixel/voxel value.Format: -cv oldValue1xoldValue2,newValue1xnewValue2.Can be used for multiple number of value changes.Defaults to 3,4."
  dicom2Nifti:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -d2n
    doc: If path to reference is present, then image comparison is done.Use '-i' to pass input DICOM image.Use '-o' to pass output image file.
hints:
  SoftwareRequirement:
    packages:
      Utilities:
        version:
          - 1.7.3.nonRelease.20190819