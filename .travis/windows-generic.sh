#!/bin/bash

# This is a Travis CI script intended to be ran in window's git bash shell

mkdir bin
cd bin

cmake -DCMAKE_INSTALL_PREFIX="./install" -DBUILD_DOCUMENTATION=OFF ..
make