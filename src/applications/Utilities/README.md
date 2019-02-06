# Utilities

Though some functionality does exist on the graphical interface, this executable is primarily designed to be used from the command line. 

## Functionality offered

1. Resizing: Done as a percentage of the original image with options for different interpolations
2. Resample: Change the input image voxel spacing
3. SanityCheck: Sanity check with a reference image to check physical characteristics
4. Information: Get information in the input image
5. Casting: Change input image type
6. UniqueValues: Get the unique values in input image
7. TestComparison: Baseline image to compare input image with
8. BoundingBox: Extracts the smallest bounding box around the mask file (optionally, can be isotropic)
9. CreateMask: Create a binary mask out of a provided (float) thresholds
10. ChangeValue: Change the specified pixel/voxel value

## To Do

- Extend image information with S/Q Form matrix information and Direction cosine of image

For new features, please put a pull or feature request via github.com/CBICA/CaPTk