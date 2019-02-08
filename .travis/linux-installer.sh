#!/bin/bash

# Download the CaPTk installer and try to do a basic install

# Any subsequent(*) commands which fail will 
# cause the shell script to exit immediately
set -e

# Binary name
BIN_NAME="CaPTk_1.6.1.bin"

# Download the latest CaPTk release
# NOTE: Please update this as the verions change!
wget --tries=inf ftp://www.nitrc.org/home/groups/captk/downloads/${BIN_NAME} --output ${BIN_NAME}

# Spawn an instance of the installer
spawn ./${BIN_NAME}

# Fake a license acceptance
expect "#? "
send "1"

# Should now install, if problems arise, travis will catch them.