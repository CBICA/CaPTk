# Scripts

This directory contains scripts relevant to CaPTk configuration, compilation, or packaging.

- `captk-pkg` - Packaging suite. Builds, installs, and packages CaPTk into an AppImage if the superbuild has already been completed.
- `linux-build` - Stripped down version of `captk-pkg` that removes the packaging step. Still implies a complete superbuild.
- `linux-build-collab` - Version of the build script meant to be used by first-time collaborators with their own external ITK, VTK, and OpenCV dependencies
- `linux-cmake-conf` - Configuration script that makes cmake install in a set of directories meant for the AppImage. Used by `captk-pkg` and thus it's derivatives.
- `linux-makeself` - The script to be ran as the installation script in the AppImage. This is for developer use only, and serves no practical purpose without the AppImage context.