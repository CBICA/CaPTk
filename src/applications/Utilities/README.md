# Utilities

Though some functionality does exist on the graphical interface, this executable is primarily designed to be used from the command line. 

## Functionality offered

- Resizing: Done as a percentage of the original image with options for different interpolations
- Resample: Change the input image voxel spacing
- SanityCheck: Sanity check with a reference image to check physical characteristics
- Information: Get information in the input image
- Casting: Change input image type
- UniqueValues: Get the unique values in input image
- TestComparison: Baseline image to compare input image with
- BoundingBox: Extracts the smallest bounding box around the mask file (optionally, can be isotropic)
- CreateMask: Create a binary mask out of a provided (float) thresholds
- ChangeValue: Change the specified pixel/voxel value
- DICOM to NIfTI
- NIfTI to DICOM
- NIfTI to DICOM-Segmentation: leverages [DCMQI](https://github.com/QIICR/dcmqi)
- Re-orient image
- Thresholding: above, below, Otsu, binary
- Coordinate conversion: from world to image and vice-versa
- Extract Image Series: convert N-D image to (N-1)-D image; useful for extract volumes from time series images
- Join Image Series: convert N-D image to (N+1)-D image; useful for constructing time series images
- Label Similarity measures: calculate various similarity statistics between two label sets
- BraTS Label Similarity measures: calculate similarity statistics between two brain tumor masks

## To Do

- Extend image information with S/Q Form matrix information and Direction cosine of image

For new features, please put a pull or feature request via github.com/CBICA/CaPTk