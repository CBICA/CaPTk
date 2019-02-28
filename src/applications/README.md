# CaPTk Applications

## Adding new C++ Applications

1. Create a folder with your application's name under `CaPTk/src/applications` *OR* instantiate a git sub-module referring to your application's repo in the same location.

2. In the root `CMakeLists.txt`, ensure you add the following lines to ensure your source files are added into CaPTk's include directory structure:

```cmake
SET( CACHED_INCLUDE_DIRS 
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/includes # change this and add values accordingly
  CACHE STRING "All include directories"
)  
```

3. 