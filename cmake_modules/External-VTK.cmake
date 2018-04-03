MESSAGE( "External project - VTK" )

SET( VTK_DEPENDENCIES )

ExternalProject_Add( VTK
  DEPENDS ${ITK_DEPENDENCIES}
  GIT_REPOSITORY ${git_protocol}://vtk.org/VTK.git
  GIT_TAG v5.10.1
  SOURCE_DIR VTK
  BINARY_DIR VTK-build
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${ep_common_args}
    -DBUILD_EXAMPLES:BOOL=OFF # examples are not needed
    -DBUILD_SHARED_LIBS:BOOL=ON # no static builds
    -DBUILD_TESTING:BOOL=OFF # testing the ITK build is not required
	-DVTK_USE_QT:BOOL=TRUE # [QT] dependency, enables better GUI
    -DVTK_USE_GUISUPPORT:BOOL=TRUE 
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/INSTALL
)

set( VTK_DIR ${CMAKE_BINARY_DIR}/VTK-build )