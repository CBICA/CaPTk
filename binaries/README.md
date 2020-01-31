# Binaries

Contains the binaries used by CI/the superbuild. These are automatically handled for you by cmake, no need to touch them.

## License Agreement

By downloading these binaries, you are assuming acceptance of the specific licenses being described in ../licenses/Combined.txt

## Qt Options

Extraction has been done using the relevant Qt installer from the web by enabling the following options for the chosen compiler:

- Charts
- Data Visualization
- WebEngine
- Designer & Developer Tools: developer-specific

## Updating The Executable of an Existing Application

1. Update the executable for one/all the platforms in their respective compressed files under ```precompiledApps```.
2. Ensure that each app has its own folder, this makes it easier for CMake to pick it up for installation, which happens in https://github.com/CBICA/CaPTk/blob/master/src/applications/CMakeLists.txt.


## Adding a New Executable

1. Add the executable for one/all the platforms in their respective compressed files under ```precompiledApps```.
2. Ensure the application has its own folder, this makes it easier for CMake to pick it up for installation, which happens in https://github.com/CBICA/CaPTk/blob/master/src/applications/CMakeLists.txt.
3. Each application can (and has) its own installation mechanism (see LIBRA and ITK-SNAP as examples)
4. The Installation control is done using **BUILD_APPNAME** CMake variable in https://github.com/CBICA/CaPTk/blob/master/CMakeLists.txt#L337, which can be toggled as needed
5. The actual installation is based on if the previous variable is enabled **AND** if the appropriate folder is found (as an example, see [LIBRA's installation](https://github.com/CBICA/CaPTk/blob/master/src/applications/CMakeLists.txt#L335))
6. Different platforms can (and do) have different installation mechanisms:
   1. Windows & Linux follow the same common path.
   2. macOS, on the other hand, needs some "special sauce" to work. 
   3. Example is DeepMedic's installation: https://github.com/CBICA/CaPTk/blob/master/src/applications/ eLists.txt#L302.
   4. LIBRA and all associated applications (Breast Texture Pipeline) do not exist on macOS.
   5. ITK-SNAP needs to be installed in a "special" way for macOS since it is distributed as an AppBundle: https:// ub.com/CBICA/CaPTk/blob/2d779d1674013cef9c7d1b1b708700461dee9a73/CMakeLists.txt#L1068.