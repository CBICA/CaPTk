# CBICA Toolkit

This project houses the common classes and functions that are used throughout CaPTk and other C++ projects under the CBICA umbrella (including collaborators). All functions/classes are cross-platform.

## Dependencies

- ITK

Changes to dependencies should be discussed in the group so that any downstream effects can be rectified.

## Functionality

- Command line parsing
- Safe Image I/O (NIfTI, DICOM)
- Image Utilities and wrappers to ITK filters
..- CreateMaskIndeces
..- GetPixelValuesFromIndeces
..- Preprocessing: histogram matching, smoothing, image comparison, orientation fix, skull stripping, resize, resample
..- Image sanity checking
..- Distance calculations: image and world coordinates
..- Create new image based on existing image
- File system functionalities
