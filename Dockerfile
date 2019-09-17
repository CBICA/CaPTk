FROM ubuntu:16.04

LABEL authors="CBICA_UPenn <software@cbica.upenn.edu>"

#general dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    mesa-utils \
    freeglut3-dev \
    wget \
    unzip \
    doxygen \
    -qq \
    gcc-5 \
    g++-5 \
    make \
    libgl-dev \
    python3-pip \
    python-numpy \
    dos2unix \
    libxkbcommon-x11-0 \
    libnss3 \
    libxcomposite-dev \
    libxcursor-dev \
    libxrender-dev \
    libxtst-dev \
    libasound2 \
    libdbus-1-dev \
    libegl1-mesa \
    libglib2.0-dev \
    libxext6 \
    libfreetype6 \
    libfreetype6-dev \
    nodejs \
    libxft-dev \
    npm; 

# installing LFS
RUN apt-get install -y sudo curl git && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
    sudo apt-get install git-lfs; \
    git lfs install;    

# installing CMake
RUN wget https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.sh; \
    mkdir /opt/cmake; \
    sh cmake-3.12.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license; \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake

# setting up the build environment
ARG GIT_LFS_SKIP_SMUDGE=1
ARG PKG_FAST_MODE=1
ARG PKG_COPY_QT_LIBS=1
ENV GIT_LFS_SKIP_SMUDGE=$GIT_LFS_SKIP_SMUDGE
ENV PKG_FAST_MODE=$PKG_FAST_MODE
ENV PKG_COPY_QT_LIBS=$PKG_COPY_QT_LIBS

# cloning CaPTk
RUN if [ ! -d "`pwd`/CaPTk" ] ; then git clone "https://github.com/CBICA/CaPTk.git"; fi 
RUN cd CaPTk &&  git pull; \
    git submodule update --init && mkdir bin

RUN cd CaPTk/bin && echo "=== Starting CaPTk Superbuild ===" && \
    if [ ! -d "`pwd`/externalApps" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/precompiledApps/linux.zip -O binaries_linux.zip; fi ; \
    if [ ! -d "`pwd`/qt" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/qt_5.12.1/linux.zip -O qt.zip; fi ; \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install_libs -Wno-dev ..

RUN cd CaPTk/bin && echo "=== Building CaPTk ===" && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install/appdir/usr/bin -Wno-dev .. && \
    make install/strip -j2;

RUN cd CaPTk && ./scripts/captk-pkg

# cleanup
RUN rm -rf CaPTk/bin/binaries_linux.zip && rm -rf CaPTk/bin/qt.zip

# define entry point
ENTRYPOINT ["/CaPTk/bin/install/appdir/usr/bin/CaPTk"]
