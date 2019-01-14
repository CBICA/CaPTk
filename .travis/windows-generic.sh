#!/bin/bash

# This is a Travis CI script intended to be ran in window's git bash shell

# Nuclear option
# rm -rf binaries
rm -rf data
rm -rf history
rm -rf src/applications/individualApps/libra/MCRInstaller.zip

mkdir bin

mv ./binaries/qt5.11.2_windows.zip ./bin/qt.zip

cd bin

cmake -DCMAKE_INSTALL_PREFIX="./install" -DBUILD_DOCUMENTATION=OFF ..
C:/Program Files (x86)/Microsoft Visual Studio/2017/Enterprise(community)/MSBuild/bin/msbuild.exe CaPTk_Depends.sln