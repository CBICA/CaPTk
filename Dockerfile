FROM cbica/captk_centos7:devtoolset-4_superbuild

LABEL authors="CBICA_UPenn <software@cbica.upenn.edu>"

RUN yum update -y

RUN yum install git

RUN cd CaPTk; \
    git pull origin master

RUN cd CaPTk/bin; \
    if [ ! -d "`pwd`/externalApps" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/precompiledApps/linux.zip -O binaries_linux.zip; fi ; \
    cmake -DITK_DIR=./bin/ITK-build -DDCMTK_DIR=./bin/DCMTK-build -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DBUILD_TESTING=OFF ..; \
    make -j2 && make install/strip; \
    cd .. && ./scripts/captk-pkg

# cleanup
RUN rm -rf CaPTk/bin/binaries_linux.zip

# define entry point
ENTRYPOINT ["/CaPTk/bin/install/appdir/usr/bin/CaPTk"]