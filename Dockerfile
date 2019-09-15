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

RUN if [ ! -d "`pwd`/cmake-3.12.4-Linux-x86_64" ] ; then wget https://cmake.org/files/v3.12/cmake-3.12.4-Linux-x86_64.tar.gz ; fi 
RUN if [ ! -d "`pwd`/cmake-3.12.4-Linux-x86_64" ] ; then tar -xvf cmake-3.12.4-Linux-x86_64.tar.gz ; fi 
RUN if [ ! -d "`pwd`/cmake-3.12.4-Linux-x86_64" ] ; then rm -rf cmake-3.12.4-Linux-x86_64.tar.gz; fi 

RUN apt-get install -y sudo curl git && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
    sudo apt-get install git-lfs; \
    git lfs install;    

ENV PATH `pwd`/cmake-3.12.4-Linux-x86_64/bin:$PATH
ENV GIT_LFS_SKIP_SMUDGE=1
ENV PKG_FAST_MODE=1
ENV PKG_COPY_QT_LIBS=1

RUN if [ ! -d "`pwd`/CaPTk" ] ; then git clone "https://github.com/CBICA/CaPTk.git"; fi \
    cd CaPTk &&  git pull; \
    rm -rf *.bin; \
    git submodule init; \
    git submodule update;

RUN cd CaPTk && echo "=== Starting CaPTk Superbuild ===" && \
    mkdir bin && cd bin && \
    if [ ! -d "`pwd`/externalApps" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/precompiledApps/linux.zip -O binaries_linux.zip; fi \
    if [ ! -d "`pwd`/qt" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/qt_5.12.1/linux.zip -O qt.zip; fi

RUN cd CaPTk/bin && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install_libs -Wno-dev ..

RUN cd CaPTk/bin && make

RUN cd CaPTk/bin && if [ ! -f "qt.zip" ] ; then rm -rf qt.zip; fi

RUN cd CaPTk/bin && \
    echo "=== Building CaPTk ===" && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install/appdir/usr/bin -Wno-dev .. && \
    make install/strip -j2;

RUN cd CaPTk && ./scripts/captk-pkg

#general dependencies
RUN apt-get update && \
    apt-get install -y sudo curl git && \
    curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
    sudo apt-get install git-lfs; \
    git lfs install; \
    if [ ! -d "`pwd`/CaPTk" ] ; then git clone "https://github.com/CBICA/CaPTk.git"; fi \
    cd CaPTk &&  git pull; \
    count=`ls -1 *.flac 2>/dev/null | wc -l`; \
    if [ $count != 0 ] ; then rm -rf *.bin; fi \
    git submodule init; \
    git submodule update; \
    echo "=== Starting CaPTk Superbuild ===" && \
    mkdir bin && cd bin && \
    if [ ! -d "`pwd`/externalApps" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/precompiledApps/linux.zip -O binaries_linux.zip; fi \
    if [ ! -d "`pwd`/qt" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/qt_5.12.1/linux.zip -O qt.zip; fi \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install_libs -Wno-dev .. && \
    make && \
    if [ ! -f "qt.zip" ] ; then rm -rf qt.zip; fi && \
    echo "=== Building CaPTk ===" && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install/appdir/usr/bin -Wno-dev .. && \
    make install/strip -j2; \
    cd .. && ./scripts/captk-pkg

# define entry point
ENTRYPOINT ["/CaPTk/bin/install/appdir/usr/bin/CaPTk"]
