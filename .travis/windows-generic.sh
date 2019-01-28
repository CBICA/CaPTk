#!/bin/bash

# This is a Travis CI script intended to be ran in window's git bash shell

# Nuclear option
# rm -rf binaries
rm -rf data
rm -rf history
rm -rf src/applications/individualApps/libra/MCRInstaller.zip

mv ./binaries/qt5.11.2_windows.zip ./qt.zip

cmake -DCMAKE_INSTALL_PREFIX="./install" -DBUILD_DOCUMENTATION=OFF .
cmake --build . --target ALL_BUILD --config Release

cmake -DCMAKE_INSTALL_PREFIX="./install" -DBUILD_DOCUMENTATION=OFF .
cmake --build . --target ALL_BUILD --config Release
