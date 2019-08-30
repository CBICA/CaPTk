FIND_PACKAGE( Qt5Core ) # not keeping version requirement here because of centos

# buld Qt from source for trusty; otherwise use pre-built binaries
IF( "${RELEASE_CODENAME}" STREQUAL "trusty" )

  OPTION( QT_BUILD_FROM_SOURCE "Build Qt ${QT_VERSION} from " OFF )
  
  FIND_PACKAGE( Qt5Core ${QT_VERSION} )
  
  IF( (NOT Qt5Core_FOUND) OR QT_BUILD_FROM_SOURCE )

    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/qtifwsilent.qs
      ${CMAKE_CURRENT_BINARY_DIR}/qtifwsilent.qs
      @ONLY 
    )
     
    EXECUTE_PROCESS(
      COMMAND curl -L -O 'https://download.qt.io/official_releases/qt/5.12/5.12.1/qt-opensource-linux-x64-${QT_VERSION}.run' && chmod +x qt-opensource-linux-x64-${QT_VERSION}.run && QT_INSTALL_DIR=/usr/local/Qt ./qt-opensource-linux-x64-${QT_VERSION}.run --platform minimal --script qtifwsilent.qs
      OUTPUT_VARIABLE TEMP
    )
    
    EXECUTE_PROCESS(
      COMMAND export PATH="/usr/local/Qt/${QT_VERSION}/gcc_64/bin/:${PATH}"
      OUTPUT_VARIABLE TEMP
    )

    EXECUTE_PROCESS(
      COMMAND export LD_LIBRARY_PATH="/usr/local/Qt/${QT_VERSION}/gcc_64/lib:${LD_LIBRARY_PATH}"
      OUTPUT_VARIABLE TEMP
    )

    EXECUTE_PROCESS(
      COMMAND export CMAKE_PREFIX_PATH="/usr/local/Qt/${QT_VERSION}/gcc_64/lib/cmake/Qt5:${CMAKE_PREFIX_PATH}"
      OUTPUT_VARIABLE TEMP
    )
    
  ENDIF()

ELSE()
  OPTION( QT_DOWNLOAD_FORCE "Force Qt binary download regardless of whether Qt was found in host machine or not" OFF )
ENDIF()

IF( (NOT Qt5Core_FOUND) OR QT_DOWNLOAD_FORCE )

  SET( FILENAME_TO_EXTRACT "qt" )
  SET( FILE_TO_EXTRACT "${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip" )
  SET( QT_EXTRACTED_DIR "${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}" )

  SET( PLATFORM_STRING "" )

  IF(WIN32)
    SET( PLATFORM_STRING "windows" )
  ELSEIF(APPLE)
    SET( PLATFORM_STRING "macos" )
  ELSE()
    SET( PLATFORM_STRING "linux" )
  ENDIF()
  
  SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/qt/${QT_VERSION}/${PLATFORM_STRING}.zip" )
  # SET( DOWNLOAD_LINK "https://github.com/CBICA/CaPTk/raw/master/binaries/qt_${QT_VERSION}/${PLATFORM_STRING}.zip" ) # the centos qt is not here, yet
  SET( LFS_FILE_TO_CHECK "${PROJECT_SOURCE_DIR}/binaries/qt_${QT_VERSION}/${PLATFORM_STRING}.zip" )

  # if(UNIX)
  #   if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  #     if(EXISTS "/etc/issue")
  #       file(READ "/etc/issue" LINUX_ISSUE)
  #       if(LINUX_ISSUE MATCHES "CentOS")
  #         SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/qt/${QT_VERSION}/${PLATFORM_STRING}.zip" )
  #       endif()
  #     endif()
  #   endif()
  # endif()

  IF( NOT EXISTS "${QT_EXTRACTED_DIR}" )
    FILE(MAKE_DIRECTORY "${QT_EXTRACTED_DIR}" )
  ENDIF()

  IF( NOT EXISTS "${FILE_TO_EXTRACT}" )
    
    # copy from LFS folder
    IF( EXISTS ${LFS_FILE_TO_CHECK} )
      CONFIGURE_FILE( ${LFS_FILE_TO_CHECK} ${FILE_TO_EXTRACT} )
    ENDIF()

    # do not re-download if the LFS fetch worked
    IF(NOT EXISTS ${FILE_TO_EXTRACT})
      MESSAGE( STATUS "Downloading pre-compiled Qt with open source license (see Qt site for more details)" )
      FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
      IF(NOT STATUS_CODE EQUAL 0)
        MESSAGE(FATAL_ERROR "Failed to download pre-compiled packages. Status=${STATUS_CODE}")
      ENDIF()
    ENDIF()

  ENDIF()

  IF( EXISTS "${FILE_TO_EXTRACT}" )

    IF( NOT EXISTS "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5" )
      MESSAGE( STATUS "Extracting pre-compiled Qt binaries" )
      EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
        WORKING_DIRECTORY ${QT_EXTRACTED_DIR}
        RESULT_VARIABLE RESULT_CODE
      )

      IF(NOT RESULT_CODE EQUAL 0)
        MESSAGE( WARNING "Extracting the pre-compiled Qt binaries failed" )
      ENDIF()
    ENDIF()

    LIST(APPEND CMAKE_PREFIX_PATH "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5")
    LIST(APPEND CMAKE_PROGRAM_PATH  "${QT_EXTRACTED_DIR}/${QT_VERSION}/bin")
        
    IF (APPLE) # THIS QT FOLDER STRUCTURE WAS NOT REQUIRED FOR 5.11.2 
      SET( CAPTK_QT5_DIR "${QT_EXTRACTED_DIR}/${QT_VERSION}/clang_64/lib/cmake/Qt5" CACHE STRING "Qt5_DIR Path for use other builds" FORCE )
      SET( ENV{PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/clang_64/lib/cmake/Qt5/:${QT_EXTRACTED_DIR}/${QT_VERSION}/clang_64/bin:$ENV{PATH}" CACHE PATH "" FORCE)
      SET( ENV{CMAKE_PREFIX_PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/clang_64/lib/cmake/Qt5/:${QT_EXTRACTED_DIR}/${QT_VERSION}/clang_64/bin:${PROJECT_BINARY_DIR}/ITK-build/:${PROJECT_BINARY_DIR}/OpenCV-build/:${PROJECT_BINARY_DIR}/VTK-build:$ENV{CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE )
      SET( ENV{CMAKE_PROGRAM_PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/clang_64/lib/cmake/Qt5/:${QT_EXTRACTED_DIR}/${QT_VERSION}/clang_64/bin:$ENV{CMAKE_PROGRAM_PATH}" CACHE PATH "" FORCE )
    ELSE()
      SET( CAPTK_QT5_DIR "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5" CACHE STRING "Qt5_DIR Path for use other builds" FORCE )
      SET(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5/" )

      SET( ENV{CMAKE_PREFIX_PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5/:${QT_EXTRACTED_DIR}/${QT_VERSION}/bin:${PROJECT_BINARY_DIR}/ITK-build/:${PROJECT_BINARY_DIR}/OpenCV-build/:${PROJECT_BINARY_DIR}/VTK-build:$ENV{CMAKE_PREFIX_PATH}" CACHE PATH "" FORCE )
      SET( ENV{LD_LIBRARY_PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/:$ENV{LD_LIBRARY_PATH}" CACHE PATH "" FORCE )
      SET( ENV{CMAKE_PROGRAM_PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5/:${QT_EXTRACTED_DIR}/${QT_VERSION}/bin:$ENV{CMAKE_PROGRAM_PATH}" CACHE PATH "" FORCE )
      SET( ENV{PATH} "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5/:${PROJECT_BINARY_DIR}/qt/${QT_VERSION}/bin:$ENV{PATH}" CACHE PATH "" FORCE)
    ENDIF()

    MESSAGE( AUTHOR_WARNING "=== [DEBUG] ENV{CMAKE_PREFIX_PATH}=$ENV{CMAKE_PREFIX_PATH}" )
    # SET(ENV{CMAKE_PREFIX_PATH} "${CMAKE_PREFIX_PATH};${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5/" )
    # SET(ENV{CMAKE_PROGRAM_PATH} "${CMAKE_PREFIX_PATH};${QT_EXTRACTED_DIR}/${QT_VERSION}/bin" )
    
    FIND_PACKAGE(Qt5 COMPONENTS Core Gui Svg Widgets WebView WebEngine WebEngineCore)

  ENDIF()
    
  LINK_DIRECTORIES(${QT_LIBRARY_DIR})
  SET( Qt5_DIR ${Qt5_DIR} CACHE STRING "Qt5_DIR Path for use in other builds" FORCE )
  # SET( CAPTK_QT5_DIR "${QT_EXTRACTED_DIR}/${QT_VERSION}/lib/cmake/Qt5" CACHE STRING "Qt5_DIR Path for use other builds" FORCE )
  LIST(APPEND CMAKE_PREFIX_PATH "${Qt5_DIR}")

ENDIF()
# set (qt_depends)
# set (qt_options)
# if (APPLE)
  # list (APPEND qt_options
    # -sdk ${CMAKE_OSX_SDK}
    # -qt-libpng)
# elseif (UNIX)
  # list (APPEND qt_depends freetype fontconfig png)
  # list (APPEND qt_options
    # -qt-xcb
    # -system-libpng
    # -I <INSTALL_DIR>/include/freetype2
    # -I <INSTALL_DIR>/include/fontconfig)
# endif()

# if (NOT WIN32)
  # list(APPEND qt_depends
    # zlib)
  # list(APPEND qt_options
    # -no-alsa
    # -no-pulseaudio
    # -system-zlib)
# else ()
  # list(APPEND qt_options
    # -qt-zlib)
# endif ()

# set(qt_EXTRA_CONFIGURATION_OPTIONS ""
    # CACHE STRING "Extra arguments to be passed to Qt when configuring.")

# set(qt_configure_command <SOURCE_DIR>/configure)
# if (WIN32)
  # set(qt_configure_command <SOURCE_DIR>/configure.bat)
# endif ()

# set(qt_build_commands)
# if (WIN32)
  # list(APPEND qt_build_commands
    # BUILD_COMMAND ${NMAKE_PATH}
    # INSTALL_COMMAND ${NMAKE_PATH} install)
# endif ()

# add_external_project_or_use_system(
    # qt
    # DEPENDS zlib ${qt_depends}
    # CONFIGURE_COMMAND
      # ${qt_configure_command}
      # -prefix <INSTALL_DIR>
      # -opensource
      # -release
      # -confirm-license
      # -nomake examples
      # -skip qtconnectivity
      # -skip qtlocation
      # -skip qtmultimedia
      # -skip qtquick1
      # -skip qtsensors
      # -skip qtserialport
      # -skip qtsvg
      # -skip qtwayland
      # -skip qtwebchannel
      # -skip qtwebengine
      # -skip qtwebkit
      # -skip qtwebsockets
      # -no-dbus
      # -no-openssl
      # -qt-libjpeg
      # -qt-pcre
      # -I <INSTALL_DIR>/include
      # -L <INSTALL_DIR>/lib
      # ${qt_options}
      # ${qt_EXTRA_CONFIGURATION_OPTIONS}
    # ${qt_build_commands}
# )

# add_extra_cmake_args(
  # -DPARAVIEW_QT_VERSION:STRING=5
  # -DQt5_DIR:PATH=<INSTALL_DIR>/lib/cmake/Qt5
# )