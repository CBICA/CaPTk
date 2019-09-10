FROM ubuntu:16.04

LABEL authors="CBICA_UPenn <software@cbica.upenn.edu>"

# update
RUN apt-get update  

#general dependencies
RUN apt-get install -y \
    build-essential \
    mesa-utils \
    freeglut3-dev \
    wget \
    git-core \
    unzip \
    doxygen \
    -qq \
    gcc-4.8 \
    g++-4.8 \
    make \
    libgl-dev \
    python3-pip \
    python-numpy \
    dos2unix \
    libxkbcommon-x11-0 \
    nodejs \
    npm 

# testing git installation
RUN add-apt-repository ppa:git-core/ppa; \
    apt-get update; \
    apt-get install -y git

# git lfs
RUN apt-get update && \
    apt-get install -y sudo curl git && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
    sudo apt-get install git-lfs; \
    git lfs install

# cmake installation
RUN wget https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz; \
    tar -xf cmake-3.12.4-Linux-x86_64.tar.gz 

# setup environment
ENV PATH=`pwd`/cmake-3.14.3-Linux-x86_64/bin \
    GIT_LFS_SKIP_SMUDGE=1 

# clone the current repo
RUN git clone https://github.com/CBICA/CaPTk.git

RUN which cmake

# start superbuild and then build CaPTk
RUN cd CaPTk && \
    echo "=== Starting CaPTk Superbuild ===" && \
    mkdir bin && cd bin && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install_libs -Wno-dev .. && \
    make -j4 && \
    echo "=== Building CaPTk ===" && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install -Wno-dev .. && \
    make install/strip -j4

# define entry point
ENTRYPOINT ["/CaPTk/bin/install/bin/CaPTk"]
