#!/bin/bash

# Download the CaPTk installer and try to do a basic install

# Any subsequent(*) commands which fail will 
# cause the shell script to exit immediately
set -e

# Download the latest CaPTk release
# NOTE: Please update this as the verions change!
wget ftp://www.nitrc.org/home/groups/captk/downloads/CaPTk_1.6.1.bin

# Spawn an instance of the installer
spawn ./CaPTk_1.6.1.bin

# Fake a license acceptance
expect "#? "
send "1"

# Should now install, if problems arise, travis will catch them.