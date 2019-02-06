ExternalProject_Add( 
  OpenCV_Contrib
  URL https://github.com/opencv/opencv_contrib/archive/3.4.5.zip
  SOURCE_DIR OpenCV_Contrib
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
)

SET( OPENCV_CONTRIB_PATH ${CMAKE_BINARY_DIR}/OpenCV_Contrib CACHE STRING "OpenCV Contrib Dir" )