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
rm -rf data
rm -rf history
rm -rf src/applications/individualApps/libra/MCRInstaller.zip

# Create binary directory
echo "[:] Creating binary directory..."
mkdir -p bin

# Move OS specific qt lib in
mv ./binaries/qt5.11.2_linux.zip ./bin/qt.zip

# Move externalApps into bin to trick CMake
mv ./binaries/externalApps.zip ./bin/

# Remove all other blobs
rm -rf binaries

# FTP override
# wget -t inf ftp://www.nitrc.org/home/groups/captk/downloads/qt/5.11.2/linux.zip
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

# # Cmake
# echo "[:] Running cmake command..."
# bash ../scripts/linux-cmake-conf

# # Make
# echo "[:] Building CaPTk..."
# make -j2

echo "[:] Done. Project has been built"