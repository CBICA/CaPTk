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
- Generic functionalities: file/directory check, sym-link creation, sub-directories, filename splitting, create/delete/copy directory, file/folder size, normalize path, get current working and/or executable directory, get/set environment variables, CSV parsing, confusion matrix, ROC, randn, 
- First order statistics using a standard vector

## Usage 

### Root CMakeLists File:

```cmake
# for projects within CaPTk/src/applications
ADD_SUBDIRECTORY( ${PROJECT_SOURCE_DIR}/../../cbica_toolkit ) 

# for projects that add CaPTk as a sub-module
ADD_SUBDIRECTORY( ${CaPTk_base_sources}/src/cbica_toolkit ) 

INCLUDE_DIRECTORIES( ${CACHED_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} CaPTk_CBICATK )
```

### Source code:

```cpp
#include "cbicaUtilities.h" // or the relevant header file

int main()
{
  // use functions as needed
}
```