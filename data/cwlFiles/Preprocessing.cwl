cwlVersion: v1.0
class: CommandLineTool
baseCommand: Preprocessing
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
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Image for processing.
  maskImage:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -m
    doc: Input Mask for processing.
  outputImage:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -o
    doc: Output Image for processing.
  histoMatch:
    type: File?
    label: NIfTI Target
    inputBinding:
      position: 1
      prefix: -hi
    doc: Match inputImage with the file provided in with this parameter.Pass the target image after '-hi'.
  hMatchBins:
    type: int?
    label: 1-1000
    inputBinding:
      position: 1
      prefix: -hb
    doc: Number of histogram bins for histogram matching.Only used for histoMatching.Defaults to 100.
  hMatchQnts:
    type: int?
    label: 1-1000
    inputBinding:
      position: 1
      prefix: -hq
    doc: Number of quantile values to match for histogram matching.Only used for histoMatching.Defaults to 40.
  zScoreNorm:
    type: boolean?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -zn
    doc: Z-Score normalization.
  zNormQuant:
    type: float?
    label: 0-100
    inputBinding:
      position: 1
      prefix: -zq
    doc: "The Lower-Upper Quantile range to remove.Default: 5,95."
  zNormCut:
    type: float?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -zc
    doc: "The Lower-Upper Cut-off (multiple of stdDev) to remove.Default: 3,3."
  n3BiasCorr:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n3
    doc: "Runs the N3 bias correction.Optional parameters: mask or bins, iterations, fitting levels."
  n3BiasIter:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n3I
    doc: Number of iterations the algorithm needs to run for.Defaults to 50.
  n3BiasFitLevl:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n3F
    doc: Number of fitting levels the algorithm needs obtain.Defaults to 4.
  n3BiasBins:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n3B
    doc: If no mask is specified, N3 Bias correction makes one using Otsu.This parameter specifies the number of histogram bins for Otsu.Defaults to 200.
  n3BiasWidth:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n3W
    doc: The full width at half maximum.Defaults to 0.150000.
  n4BiasCorr:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n4
    doc: "Runs the N4 bias correction.Optional parameters: mask or bins, iterations, fitting levels."
  n4BiasIter:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n4I
    doc: Number of iterations the algorithm needs to run for.Defaults to 50.
  n4BiasFitLevl:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n4F
    doc: Number of fitting levels the algorithm needs obtain.Defaults to 4.
  n4BiasBins:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n4B
    doc: If no mask is specified, N3 Bias correction makes one using Otsu.This parameter specifies the number of histogram bins for Otsu.Defaults to 200.
  n4BiasWidth:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n4W
    doc: The full width at half maximum.Defaults to 0.150000.
  susanSmooth:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ss
    doc: Susan smoothing of an image.
  susanSigma:
    type: float?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ssS
    doc: Susan smoothing Sigma.Defaults to 0.500000.
  susanRadius:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ssR
    doc: Susan smoothing Radius.Defaults to 1.
  susanThresh:
    type: float?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ssT
    doc: Susan smoothing Intensity Variation Threshold.Defaults to 80.000000.
  p1p2norm:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -p12
    doc: P1-P2 normalization required for skull stripping.
  debugMode:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -d
    doc: "Enabled debug mode.Default: 0."
hints:
  SoftwareRequirement:
    packages:
      Preprocessing:
        version:
          - 1.7.3.nonRelease.20190819