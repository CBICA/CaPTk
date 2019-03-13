# CaPTk Applications

## Adding new C++ Applications

1. Create a folder with your application's name under `CaPTk/src/applications` *OR* instantiate a git sub-module referring to your application's repo in the same location.

2. In the root `CMakeLists.txt`, ensure you add the following lines to ensure your source files are added into CaPTk's include directory structure:

```cmake
CMAKE_MINIMUM_REQUIRED(VERSION 3.7.2)

# this is needed to ensure that project addition is done in accordance with how CaPTk expects
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake_modules/)
INCLUDE( CaPTk_macros )

CAPTK_ADD_PROJECT( ${YOUR_AWESOME_PROJECT_NAME} ${YOUR_AWESOME_PROJECT_VERSION} )

##
# set your project up here
##

# add the executable using macro so that packaging is done properly 
CAPTK_ADD_EXECUTABLE( ${YOUR_AWESOME_PROJECT_NAME} 
  ${PROJECT_SOURCE_DIR}/src/yourAwesomeProject.cxx 
  ${LIBNAME_CBICATK} ${ITK_LIBRARIES} ${OpenCV_LIBRARIES} )

# ensure all include directories (including nested dependencies) are captured for CaPTk 
SET( CACHED_INCLUDE_DIRS 
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/includes # change this and add values accordingly
  CACHE STRING "All include directories"
)  
```

3. 