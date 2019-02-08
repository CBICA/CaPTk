#!/bin/bash

# Download the CaPTk installer and try to do a basic install

# Any subsequent(*) commands which fail will 
# cause the shell script to exit immediately
set -e

# Binary name
BIN_NAME="CaPTk_1.6.1.bin"

# Download the latest CaPTk release
# NOTE: Please update this as the verions change!
# wget --tries=inf ftp://www.nitrc.org/home/groups/captk/downloads/${BIN_NAME}

# Make bin executable
chmod +x ${BIN_NAME}

# Spawn an instance of the installer and fake an acceptance
expect <<END
    spawn ./${BIN_NAME}
    expect "#? \r"
    send -- "1\r"
    expect eof
END

# Should now install, if problems arise, travis will catch them.