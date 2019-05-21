FROM ubuntu:16.04

MAINTAINER CBICA_UPenn software@cbica.upenn.edu

#update
RUN apt-get update && \
    apt-get install -y sudo curl git
    #curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
    #sudo apt-get install git-lfs

#general dependencies
RUN apt-get install -y \
    build-essential \
    mesa-common-dev \
    freeglut3-dev \
    wget \
    git-core \
    unzip \
    doxygen \
    gcc \
    g++ \
    make \
    libgl-dev \
    python3-pip \
    python-numpy \
    dos2unix \
    libxkbcommon-x11-0 \
    libxt-dev \
    libglib2.0-0
    
#RUN git lfs install

RUN ln -s `locate libc.so.6` /lib/libc.so

# install latest cmake 
RUN wget https://github.com/Kitware/CMake/releases/download/v3.14.3/cmake-3.14.3-Linux-x86_64.tar.gz && \
    tar -xzf cmake-3.14.3-Linux-x86_64.tar.gz

# clone the current repo
RUN git clone https://github.com/CBICA/CaPTk.git

# start superbuild and then build CaPTk
RUN export PATH=/cmake-3.14.3-Linux-x86_64/bin:$PATH && \
    which cmake && \
    cd CaPTk && \
    echo "=== Starting CaPTk Superbuild ===" && \
    mkdir bin && cd bin && \
    cmake -DCMAKE_INSTALL_PREFIX=./install/appdir/usr \
    -DQT_DOWNLOAD_FORCE=ON \
    -Wno-dev .. && \
    make && \
    echo "=== Building CaPTk ===" && \
    cmake -DCMAKE_INSTALL_PREFIX=./install/appdir/usr \
    -DQT_DOWNLOAD_FORCE=ON \
    -Wno-dev .. && \
    make && \
    make install/strip

# define entry point
#ENTRYPOINT ["/CaPTk/bin/CaPTk"]
ENTRYPOINT ["/CaPTk/bin/install/bin/CaPTk"]
