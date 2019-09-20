FROM cbica/captk_centos7:devtoolset-4_superbuild

LABEL authors="CBICA_UPenn <software@cbica.upenn.edu>"

RUN yum update -y

RUN yum install git

# RUN cd CaPTk; \
#     git remote add myFork https://github.com/sarthakpati/CaPTk.git; \
#     git remote -v; \
#     git fetch myFork; \
#     git checkout -b add-other-changes; \
#     git pull myFork master; \
#     git push origin HEAD; \
#     git submodule init && git submodule update; \
#     git log -1

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