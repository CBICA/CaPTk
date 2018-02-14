set (qt_depends)
set (qt_options)
if (APPLE)
  list (APPEND qt_options
    -sdk ${CMAKE_OSX_SDK}
    -qt-libpng)
elseif (UNIX)
  list (APPEND qt_depends freetype fontconfig png)
  list (APPEND qt_options
    -qt-xcb
    -system-libpng
    -I <INSTALL_DIR>/include/freetype2
    -I <INSTALL_DIR>/include/fontconfig)
endif()

if (NOT WIN32)
  list(APPEND qt_depends
    zlib)
  list(APPEND qt_options
    -no-alsa
    -no-pulseaudio
    -system-zlib)
else ()
  list(APPEND qt_options
    -qt-zlib)
endif ()

set(qt_EXTRA_CONFIGURATION_OPTIONS ""
    CACHE STRING "Extra arguments to be passed to Qt when configuring.")

set(qt_configure_command <SOURCE_DIR>/configure)
if (WIN32)
  set(qt_configure_command <SOURCE_DIR>/configure.bat)
endif ()

set(qt_build_commands)
if (WIN32)
  list(APPEND qt_build_commands
    BUILD_COMMAND ${NMAKE_PATH}
    INSTALL_COMMAND ${NMAKE_PATH} install)
endif ()

add_external_project_or_use_system(
    qt
    DEPENDS zlib ${qt_depends}
    CONFIGURE_COMMAND
      ${qt_configure_command}
      -prefix <INSTALL_DIR>
      -opensource
      -release
      -confirm-license
      -nomake examples
      -skip qtconnectivity
      -skip qtlocation
      -skip qtmultimedia
      -skip qtquick1
      -skip qtsensors
      -skip qtserialport
      -skip qtsvg
      -skip qtwayland
      -skip qtwebchannel
      -skip qtwebengine
      -skip qtwebkit
      -skip qtwebsockets
      -no-dbus
      -no-openssl
      -qt-libjpeg
      -qt-pcre
      -I <INSTALL_DIR>/include
      -L <INSTALL_DIR>/lib
      ${qt_options}
      ${qt_EXTRA_CONFIGURATION_OPTIONS}
    ${qt_build_commands}
)

add_extra_cmake_args(
  -DPARAVIEW_QT_VERSION:STRING=5
  -DQt5_DIR:PATH=<INSTALL_DIR>/lib/cmake/Qt5
)