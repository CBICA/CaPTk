#!/bin/bash

# Any subsequent(*) commands which fail will 
# cause the shell script to exit immediately
set -e

# CaPTk Packager

# Cmake command to run from /trunk/bin
# We need this directory structure for appimages to be generated
CMAKE_CMD () {
  cmake -DCMAKE_INSTALL_PREFIX="./install" -DBUILD_DOCUMENTATION=OFF ..
}

# Process files for linux

PROCESS_CMD () {
  INSTALL_DIR=`pwd`
  # Needs to be marked as executable for LIBRA
  cd ./libexec/MCR/v84/bin/glnxa64
  chmod +x matlab_helper
  # Go back up to /bin
  cd ${INSTALL_DIR}/bin
  # dos2unix stage
  dos2unix libra
  chmod +x libra
  dos2unix ConfettiGUI.py
  chmod +x ConfettiGUI.py
  cd ../libexec
  dos2unix ConfettiCore.py
  chmod +x ConfettiCore.py
  chmod +x libra
  # Up to install
  cd ${INSTALL_DIR}
}

# linuxdeployqt command
LINDEPQT_CMD () {
  APP_DIR=`pwd`
  cd ./appdir/usr/
  cp /usr/lib/x86_64-linux-gnu/nss/* ./lib
  echo "[:] Processing files for linux..."
  PROCESS_CMD
  cd ${APP_DIR}
  # Download linuxdeployqt and run it
  echo "[:] Downloading linuxdeployqt..."
  wget https://github.com/probonopd/linuxdeployqt/releases/download/5/linuxdeployqt-5-x86_64.AppImage
  chmod +x linuxdeployqt-5-x86_64.AppImage
  ./linuxdeployqt-5-x86_64.AppImage ./appdir/usr/share/applications/*.desktop -ignore-glob=usr/{libexec/MCR/**,lib/snap-3.6.0-rc1/**} -appimage
  mv ./CaPTk*.AppImage ../../CaPTk.bin
  cd ../../
}
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
rm -rf binaries
rm -rf data
rm -rf history
rm -rf src/applications/individualApps/libra/MCRInstaller.zip

# Create binary directory
echo "[:] Creating binary directory..."
mkdir bin
cd bin

# FTP override
wget -t inf ftp://www.nitrc.org/home/groups/captk/downloads/qt/5.11.2/linux.zip
mv linux.zip qt.zip

# Cmake
echo "[:] Running cmake command to configure superbuild..."
CMAKE_CMD

ls -a

# Make install/strip
echo "[:] Building depends..."
make -j2

ls -a

cd ..

ls -a

echo "[:] Running cmake command to configure CaPTk..."
CMAKE_CMD

echo "[:] Building CaPTk..."
make -j2
