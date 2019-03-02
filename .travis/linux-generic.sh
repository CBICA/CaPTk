#!/bin/bash

# Build CaPTk on a Travis environment

# Any subsequent(*) commands which fail will 
# cause the shell script to exit immediately
set -e

# Cmake command to run from /trunk/bin
# We need this directory structure for appimages to be generated
CMAKE_CMD () {
  cmake -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DBUILD_DOCUMENTATION=OFF -DQT_DOWNLOAD_FORCE=ON ..
}

# Run the fixes for linux
FIX_CMD () {
  # Save the install dir for later
  INSTALL_DIR=`pwd`

  echo "[:] Running Linux Fixes..."

  # Fix 1: LIBRA is awful on linux, and it only gets worse from here
  # 
  # For some reason LIBRA won't get added to the bin or libexec folders
  # So just manually do it
  # In addition LIBRA needs the Matlab runtime with one VERY SPECIFIC file
  # to be marked as executable or else the whole thing WILL NOT work. We'll get to this in a second
  echo "[:] Adding LIBRA..."
  cd ${SRC_DIR}
  # cp ./src/applications/individualApps/libra/bin/libra ./bin/install/appdir/usr/bin/
  # cp ./src/applications/individualApps/libra/libexec/libra ./bin/install/appdir/usr/libexec/
  # cp -r ./src/applications/individualApps/libra/MCR ./bin/install/appdir/usr/libexec/\

  # Move into install
  cd ./bin/install/appdir/usr

  # Fix 2: WebEngine stuff
  #
  # From what I understand, linuxdeployqt will pull in all the required packages, binaries and link them dynamically
  # so that it can wrap that in an AppImage so those links are preserved at runtime.
  #
  # However, the WebEngine requires libraries that chromium uses (Since it's literally just a chromium instance)
  # but for SOME reason linuxdeployqt doesn't pull in. So we need to manually add these nss libraries in for the webengine
  # to work.
  echo "[:] Adding NSS libs..."
  # cp /usr/lib/x86_64-linux-gnu/nss/* ./lib

  # Fix 3: LIBRA was made on windows and none of the stuff is marked as executable by default
  # All of these need to be exec or else it flat out doesnt work.
  # Needs to be marked as executable for LIBRA

  echo "[:] Marking MATLAB binaries as executable..."

  # The magic file that we need from earlier
  # cd ./libexec/MCR/v84/bin/glnxa64
  # chmod +x matlab_helper

  # Go back up to /bin
  cd ${INSTALL_DIR}/bin
  
  # Fix 4: DOS sucks and that file format literally cannot be run on linux so an external tool is needed
  # to mark them as usable

  echo "[:] Running dos2unix..."

  # Libra needs it
  dos2unix libra
  chmod +x libra
  # So does confetti
  dos2unix ConfettiGUI.py
  chmod +x ConfettiGUI.py
  # And also the files in libexec
  cd ../libexec
  dos2unix ConfettiCore.py
  chmod +x ConfettiCore.py
  chmod +x libra

  # Go back up to install
  cd ${INSTALL_DIR}
}

# Linuxdeployqt packaging command
LINDEPQT_CMD () {
  # Download linuxdeployqt and run it
  echo "[:] Downloading linuxdeployqt..."
  wget --tries=inf https://github.com/probonopd/linuxdeployqt/releases/download/5/linuxdeployqt-5-x86_64.AppImage
  chmod +x linuxdeployqt-5-x86_64.AppImage
  # Okay here's where it gets insane
  # So we need to run linuxdeployqt, and show it where our *.desktop files are so it can generate an AppImage that
  # properly can provide a direct link to the CaPTk binary within the /bin folder of the appdir structure
  # 
  # After that it attempts to link EVERYTHING IT CAN FIND
  # Yes, that means at some point it will try to link the **ENTIRE** MATLAB runtime and ITK-snap libraries even though they don't need it
  # So this insane regex was created to make sure linuxdeployqt just flat out recursively ignores those directories
  ./linuxdeployqt-5-x86_64.AppImage ./appdir/usr/share/applications/*.desktop -ignore-glob=usr/{libexec/MCR/**,lib/snap-3.6.0-rc1/**} -appimage
  # Move the appimage up to the root
  mv ./CaPTk*.AppImage ../../CaPTk.bin
  # And likewise move there as well
  cd ../../
}

# Store the root dir this is being run from
SRC_DIR=`pwd`

echo "[!] Warn: This script is intended for CI use. Only use it if you know what you are doing."


echo "[:] Starting CaPTk packaging process..."

echo "[?] Checking if you are in trunk..."
# Test to see if the user is in trunk
# First see if CMakeLists.txt exists
if [[ -e CMakeLists.txt ]] ; then
  # See if it contains PROJECT( CaPTk )
  if [[ -z  `cat CMakeLists.txt | grep "PROJECT( CaPTk )"`  ]] ; then
    echo "[!] Error: You do not appear to be within trunk of CaPTk (Project is not CaPTk)"
    exit -1
  fi
else
  echo "[!] Error: You do not appear to be in trunk (CMakeLists.txt not found)"
  exit -1
fi

# Nuclear option
# rm -rf binaries
# rm -rf data
rm -rf history
rm -rf src/applications/individualApps/libra/MCRInstaller.zip

# Create binary directory
echo "[:] Creating binary directory..."
mkdir -p bin

# Move OS specific qt lib in
mv ./binaries/qt5.12.1_linux.zip ./bin/qt.zip

# Move externalApps into bin to trick CMake
mv ./binaries/externalApps.zip ./bin/

# Remove all other blobs
rm -rf binaries

# FTP override
# wget -t inf ftp://www.nitrc.org/home/groups/captk/downloads/qt/5.12.1/linux.zip
# mv linux.zip qt.zip

cd bin

# Extract externalApps
unzip externalApps.zip &> /dev/null

# Create test data dir to skip ftp download
mkdir -p testing
mkdir -p ./testing/TestData

# Cmake
echo "[:] Running cmake command to configure superbuild... [stdout omitted]"
CMAKE_CMD #&> /dev/null

# Make install/strip
echo "[:] Building depends... [stdout omitted]"
make -j2 #&> /dev/null

# Build a fast installation
cd ..

echo "[:] Starting CaPTk building process..."

echo "[?] Checking if you are in trunk..."
# Test to see if the user is in trunk
# First see if CMakeLists.txt exists
if [[ -e CMakeLists.txt ]] ; then
  # See if it contains PROJECT( CaPTk )
  if [[ -z  `cat CMakeLists.txt | grep "PROJECT( CaPTk )"`  ]] ; then
    echo "[!] Error: You do not appear to be within trunk of CaPTk (Project is not CaPTk)"
    exit -1
  fi
else
  echo "[!] Error: You do not appear to be in trunk (CMakeLists.txt not found)"
  exit -1
fi

# Enter binary directory
echo "[:] Entering binary directory..."
cd bin

# Cmake
echo "[:] Running cmake command..."
bash ../scripts/linux-cmake-conf

# Make
echo "[:] Building CaPTk..."
make -j2

ls ./src/applications/FeatureExtraction/

cp FeatureExtraction.cwl ./src/applications/FeatureExtraction/

make -j2 install/strip

# # Fix
# cd install
# APP_DIR=`pwd`
# cd ./appdir/usr/
# # We need to fix the project first
# FIX_CMD
# cd ${APP_DIR}

# # Build package
# echo "[:] Building Package..."
# LINDEPQT_CMD

# # Step 1: Get makeself, and extract it
# wget --tries=inf https://github.com/megastep/makeself/releases/download/release-2.4.0/makeself-2.4.0.run
# chmod +x ./makeself-2.4.0.run
# ./makeself-2.4.0.run
# rm makeself-2.4.0.run

# # Step 2: Create a package directory for makeself to package
# mkdir pkg
# # We package the AppImage, the license, and a script to display the agreement and install the AppImage
# mv CaPTk.bin ./pkg/
# mv ./bin/install/appdir/usr/bin/*.cwl ./pkg/
# cp ./licenses/Linux-Combined.txt ./pkg/
# cp ./scripts/linux-makeself ./pkg/
# chmod +x ./pkg/linux-makeself

# # Step 3: Wrap it all upp

# # Version number for CaPTk
# ver=`./pkg/CaPTk.bin -v | grep Version | awk '{print $2}'`

# # Create the installer
# ./makeself-2.4.0/makeself.sh --gzip pkg/ CaPTk-${ver}-Installer.bin "CaPTk Linux Installer" ./linux-makeself

# # Cleanup
# rm -rf ./pkg
# rm -rf ./makeself-2.4.0*

echo "[:] Done. Project has been built"