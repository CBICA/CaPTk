FROM cbica/captk_centos7:devtoolset-4_superbuild

LABEL authors="CBICA_UPenn <software@cbica.upenn.edu>"

RUN yum update -y

RUN yum install git

# We will do git pull on the CBICA/CaPTk master, because that is the repo using which the base image is made
# We will not do compiles on the PR because the idea is that the Xenial build will check the build status of
# the PR in any case.
RUN cd CaPTk; \ 
    git pull origin master

RUN cd CaPTk/bin; \
    if [ ! -d "`pwd`/externalApps" ] ; then wget https://github.com/CBICA/CaPTk/raw/master/binaries/precompiledApps/linux.zip -O binaries_linux.zip; fi ; \
    cmake -DITK_DIR=./bin/ITK-build -DDCMTK_DIR=./bin/DCMTK-build -DCMAKE_INSTALL_PREFIX="./install/appdir/usr" -DBUILD_TESTING=OFF ..; \
    make && make install/strip; 
    #cd .. && ./scripts/captk-pkg

# cleanup
RUN rm -rf CaPTk/bin/binaries_linux.zip

# set up the docker for GUI
ENV QT_X11_NO_MITSHM=1
ENV QT_GRAPHICSSYSTEM="native"

# define entry point
ENTRYPOINT ["/work/CaPTk/bin/install/appdir/usr/bin/CaPTk"]
