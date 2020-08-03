# Unit Testing

## About

CaPTk's Unit Testing framework is build around [CTest](https://cmake.org/cmake/help/latest/manual/ctest.1.html)

## Running Unit Tests

### Windows

Run the **RUN_TEST** project on Visual Studio

### GCC/LLVM 

Run ```make test``` command on the build folder.

## Adding Tests

### Adding a new Test that invokes an existing Executable

1. Open https://github.com/CBICA/CaPTk/blob/master/testing/CMakeLists.txt.
2. Add Test via [the standard CTest format](https://cmake.org/cmake/help/v3.0/command/add_test.html).
3. This will get invoked when the developer invokes the "running of the tests".

### Adding a Test that does Functional Testing

1. Add the class/function that needs to be tested in the file https://github.com/CBICA/CaPTk/blob/master/testing/AllTests.cxx.
2. Ensure you have the testing data (not needed if creating synthetic data).

## Updating/Adding Testing Data

- Data is present in the NITRC FTP site (same place has the installers): https://github.com/CBICA/CaPTk/blob/master/testing/CMakeLists.txt#L20.
- Perform same tasks as "updating/adding" pre-compiled binaries: https://github.com/CBICA/CaPTk/blob/master/binaries/README.md
- Send updated ZIP to the FTP site (ensure you have write privileges beforehand).

### Wish-List

- Move the ```TestData.zip``` to LFS for clarity.