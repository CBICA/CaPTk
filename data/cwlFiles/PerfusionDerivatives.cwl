cwlVersion: v1.0
class: CommandLineTool
baseCommand: PerfusionDerivatives
inputs:
  input:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -i
    doc: The input DSC-MRI image..
  echo time:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -e
    doc: The echo time in seconds..
  PSR:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -p
    doc: The Percent Signal Recovery image (1=YES, 0=NO, 1 (Default)).
  peakHeight:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -pH
    doc: The Peak Height image (1=YES, 0=NO, 1 (Default)).
  apRCBV:
    type: string
    label: none
    inputBinding:
      position: 1
      prefix: -r
    doc: Automatially-extracted proxy to reletive cerebral volume image (1=YES, 0=NO, 1 (Default)).
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
      PerfusionDerivatives:
        version:
          - 1.7.3.nonRelease.20190819