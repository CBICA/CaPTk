FROM cbica/captk_ubuntu:1604superbuild

LABEL authors="CBICA_UPenn <software@cbica.upenn.edu>"

# Future superbuilds should perform this step ahead of time so we don't have to do it here, but it's fairly lightweight
RUN apt-get install -y doxygen texlive

RUN cd /work/CaPTk && \
    git pull origin master

RUN cd /work/CaPTk/bin; \
    if [ ! -d "/work/CaPTk/bin/externalApps" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/precompiledApps/linux.zip -O binaries_linux.zip$
    cmake -DITK_DIR=./bin/ITK-build -DDCMTK_DIR=./bin/DCMTK-build -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DBUILD_TESTING=OFF ..; \
    make && make install/strip && rm -rf /work/CaPTk/bin/binaries_linux.zip;
    #cd .. && ./scripts/captk-pkg


RUN apt-get install -y libegl1-mesa
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/work/CaPTk/bin/install/appdir/usr/bin/CaPTk"]

