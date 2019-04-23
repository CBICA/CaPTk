FROM ubuntu:16.04

MAINTAINER CBICA_UPenn software@cbica.upenn.edu

#update
RUN apt-get update && \
    apt-get install -y sudo curl git && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
    sudo apt-get install git-lfs

#general dependencies
RUN apt-get install -y \
    wget \
    cmake \
    git-core \
    unzip \
    doxygen 
    
RUN git lfs install

# clone the current repo
RUN git clone https://github.com/CBICA/CaPTk.git

# start superbuild and then build CaPTk
RUN cd CaPTk; \
    echo "=== Starting CaPTk Superbuild ===" \
    mkdir bin; cd bin; \
    cmake -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=./install_libs \
    -Wno-dev ..; \
    make; \
    echo "=== Building CaPTk ===" \
    cmake -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=./install \
    -Wno-dev ..; \
    make install/strip

# define entry point
ENTRYPOINT ["/CaPTk/bin/install/bin/CaPTk"]