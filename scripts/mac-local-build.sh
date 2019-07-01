#!/bin/bash

# CaPTk Packager
# Cmake command to run from /trunk/bin
# We need this directory structure for appimages to be generated
CAPTK_CMD () {
# rm -rf *

# export CC=""
# export CXX=""
# export LDFLAGS=""
# export CPPFLAGS=""

# #git lfs install && git lfs fetch --all

# ### COMMENT OUT THE 3 LINES BELOW IF DEPENDENCY MANAGER HAS BEEN BUILT
# echo "Run Dependency Manager"
# echo $CC

# export CMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE
# export CMAKE_PREFIX_PATH="/Library/TeX/texbin"

# cmake ../ -DCMAKE_INSTALL_PREFIX="./superbuild"
# cmake ../ -DCMAKE_INSTALL_PREFIX="./superbuild"

# make -j 2

# rm CMakeCache.txt
sudo rm -rf CaPTk_*

export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++
export LDFLAGS="-L/usr/local/opt/llvm/lib"
export CPPFLAGS="-L/usr/local/opt/llvm/include"

echo "Run CaPTk Build"

cmake ../
cmake ../
make -j 2

version=$(grep -i -e "project_version:*" CMakeCache.txt | cut -c24-)
pkgname="_Installer"
pkgname="$version$pkgname"

sudo rm -rf CaPTk_*.app/Contents/Resources/bin/*.app

make -j 2

rm -rf *.pkg 
rm -rf _CPack*
make package

# chmod +x CaPTk_$version.app/Contents/Resources/bin/dcm2nii
# chmod +x CaPTk_$version.app/Contents/Resources/bin/ITK-SNAP.app/Contents/MacOS/ITK-SNAP

pkgbuild --version $version --identifier com.cbica.captk --install-location /Applications --component ./_CPack_Packages/OSX/DragNDrop/CaPTk_$version/CaPTk_$version.app/  ./CaPTk_$version.pkg

productbuild --synthesize --package CaPTk_$version.pkg ./distribution.xml

xml='<?xml version="1.0" encoding="utf-8"?>
<installer-gui-script minSpecVersion="1">
    <title>CaPTk_'"$version"'</title>
    <license file="Combined.txt"></license>
    <pkg-ref id="com.cbica.captk"/>
    <options customize="never" require-scripts="false"/>
    <choices-outline>
        <line choice="default">
            <line choice="com.cbica.captk"/>
        </line>
    </choices-outline>
    <choice id="default"/>
    <choice id="com.cbica.captk" visible="false">
        <pkg-ref id="com.cbica.captk"/>
    </choice>
    <pkg-ref id="com.cbica.captk" version="$version" onConclusion="none">CaPTk_'"$version"'.pkg</pkg-ref>
</installer-gui-script>' 

echo $xml > "./distribution.xml"

productbuild --distribution ./distribution.xml --resources ./_CPack_Packages/OSX/DragNDrop/CaPTk_$version/CaPTk_$version.app/Contents/Resources/license/ --package-path . ./CaPTk_$pkgname.pkg

}

cd ../

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

# Create binary directory
echo "[:] Creating binary directory..."
mkdir -p bin

cd bin

# Cmake
echo "[:] Running cmake command and build CaPTk..."
CAPTK_CMD

echo "[:] Done. Built test target"
