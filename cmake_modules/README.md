# Cmake_modules

Contains all of the CMake files that are included via other CMakeLists.txt. If you are looking to add individual functions or macros, don't create an additional file here, please use `CaPTk_macros.cmake` and add them there.

If you'd like to change what happens during the post install step of `make install`, add your changes to `PostInstall.cmake`

## Linuxdeployqt ?

You may have noticed that there are `Macdeployqt.cmake` and `Windeployqt.cmake` but no `Linuxdeployqt.cmake`. Qt does not have a built in version of linuxdeployqt (as of 5.11) that works for this project. We instead perform this stage via `captk-pkg` utilizing [linuxdeployqt](https://github.com/probonopd/linuxdeployqt).

## Superbuild

This folder will contain the SuperBuilds (the files in links need to be updated) for all the dependencies. The mechanism is as follows:

- Check if ITK is present; if yes, then proceed to attempt a CaPTk build.
- If ITK was not found or the CMake variable *CAPTK_SUPERBUILD_FORCE* was set, then directly invoke `Superbuild.cmake`.
- If Qt was not found in the system path or if the CMake variable *QT_DOWNLOAD_FORCE* was set, download the appropriate Qt binaries (currently only supported for Windows, Linux and macOS) using the open source Qt license.