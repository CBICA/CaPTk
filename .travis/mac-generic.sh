#!/bin/bash

# CaPTk Packager
# Cmake command to run from /trunk/bin
# We need this directory structure for appimages to be generated
CAPTK_CMD () {
export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++
export LDFLAGS="-L/usr/local/opt/llvm/lib"
export CPPFLAGS="-L/usr/local/opt/llvm/include"

export CMAKE_PREFIX_PATH=/Library/TeX/texbin
export CMAKE_PREFIX_PATH="${TRAVIS_BUILD_DIR}/bin/ITK-build:$CMAKE_PREFIX_PATH"

rm -rf /usr/local/opt/qt
rm -rf /usr/local/Cellar/qt
cp -r qt/5.11.2 /usr/local/opt/qt
cp -r qt /usr/local/Cellar/qt
export CMAKE_PREFIX_PATH="/usr/local/opt/qt/lib/cmake/Qt5:/usr/local/opt/qt/bin:$CMAKE_PREFIX_PATH"

ls /usr/local/opt/qt/

cmake ../
echo "Run Dependency Manager"
make & sleep 5800; kill $! 
# make

# echo "Run CaPTk Build"
# cmake ../
# cmake ../
# make

}

###########################
echo "[!] Warn: This script is intended for CI use. Only use it if you know what you are doing."

echo "[:] Starting CaPTk packaging process..."

echo "[?] Checking if you are in trunk..."
# Test to see if the user is in trunk
# First see if CMakeLists.txt exists
if [[ -e CMakeLists.txt ]] ; then
# See if it contains PROJECT( CaPTk )
if [[ -z `cat CMakeLists.txt | grep "PROJECT( CaPTk )"` ]] ; then
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
# mv ${TRAVIS_BUILD_DIR}/binaries/qt5.11.2_macos.zip ./bin/qt.zip

# Move externalApps into bin to trick CMake
mv ${TRAVIS_BUILD_DIR}/binaries/externalApps.zip ./bin/

# Remove all other blobs
# rm -rf ${TRAVIS_BUILD_DIR}/binaries

cd bin

# Extract externalApps
unzip externalApps.zip &> /dev/null

# Create test data dir to skip ftp download
mkdir testing
mkdir ./testing/TestData

# Cmake
echo "[:] Running cmake command and build CaPTk..."
CAPTK_CMD

# # Make install/strip
# echo "[:] Building CaPTk..."
# make

echo "[:] Done. Built test target"
