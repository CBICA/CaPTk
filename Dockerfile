FROM centos:7

LABEL authors="CBICA_UPenn <software@cbica.upenn.edu>"

#update
RUN yum -y update bash

RUN yum update -y

# taken from https://github.com/sclorg/devtoolset-container/blob/master/6-toolchain/Dockerfile
RUN yum install -y centos-release-scl-rh wget && \
    INSTALL_PKGS="devtoolset-4-gcc devtoolset-4-gcc-c++ devtoolset-4-gcc-gfortran devtoolset-4-gdb devtoolset-4-toolchain" && \
    yum install -y --setopt=tsflags=nodocs $INSTALL_PKGS && \
    rpm -V $INSTALL_PKGS && \
    yum -y clean all --enablerepo='*'

RUN rpmkeys --import file:///etc/pki/rpm-gpg/RPM-GPG-KEY-CentOS-7 && \
    yum -y install centos-release-scl

# cmake installation
RUN wget https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz; \
    tar -xf cmake-3.12.4-Linux-x86_64.tar.gz 

#general dependencies
RUN yum install -y \
    sudo \
    #devtoolset-4 \
    #devtoolset-4-gcc* \
    # gcc \
    # gcc-c++ \
    make \
    yum-utils \
    wget \
    git-core \
    lapack \
    lapack-devel \
    unzip \
    tcl \
    tcl-devel \
    tk \
    tk-devel \
    fftw \
    fftw-devel \
    mpich \
    mpich-devel \
    mesa-libGL \
    mesa-libGL-devel \
    xorg-x11-server-Xorg \
    xorg-x11-xauth \
    xorg-x11-apps \
    libXt-devel \
    tcl \
    time \
    libmpc-devel \
    mpfr-devel \
    gmp-devel \
    doxygen

# LFS install
RUN yum install -y epel-release git; \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.rpm.sh | sudo bash; \
    yum install -y git-lfs; \
    git lfs install

# nodejs is needed for azure
RUN curl https://raw.githubusercontent.com/creationix/nvm/v0.34.0/install.sh | bash; \
    #curl -sL https://rpm.nodesource.com/setup_12.x | sudo -E bash -; \
    curl -sL https://rpm.nodesource.com/setup | bash -; \
    yum install -y nodejs
    #curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.33.11/install.sh | bash; \
    #nvm install -y node; \

# setup environment
ENV HOME=/opt/app-root/src \
    PATH=/opt/app-root/src/bin:/opt/app-root/bin:/opt/rh/devtoolset-4/root/usr/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:`pwd`/cmake-3.14.3-Linux-x86_64/bin \
    GIT_LFS_SKIP_SMUDGE=1 \
    NVM_DIR="$HOME/.nvm"

# clone the current repo
RUN git clone https://github.com/CBICA/CaPTk.git

# start superbuild and then build CaPTk
RUN which cmake && \
    cd CaPTk && \
    echo "=== Starting CaPTk Superbuild ===" && \
    mkdir bin && cd bin && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install_libs -Wno-dev .. && \
    make -j4 && \
    echo "=== Building CaPTk ===" && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install -Wno-dev .. && \
    make install/strip -j4

# define entry point
ENTRYPOINT ["/CaPTk/bin/install/bin/CaPTk"]
