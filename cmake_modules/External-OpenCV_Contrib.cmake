
MESSAGE( STATUS "Adding OpenCV_Contrib-${OPENCV_VERSION} ...")

ExternalProject_Add( 
  OpenCV_Contrib
  URL https://github.com/opencv/opencv_contrib/archive/${OPENCV_VERSION}.zip
  SOURCE_DIR OpenCV_Contrib
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND cmake -E echo "Skipping install step."
)

SET( OPENCV_CONTRIB_PATH ${CMAKE_BINARY_DIR}/OpenCV_Contrib CACHE STRING "OpenCV Contrib Dir" )