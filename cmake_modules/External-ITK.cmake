## Superbuild module for building ITK externally.

MESSAGE( "External project - ITK" )

SET( ITK_DEPENDENCIES )

IF(UNIX)
  
  ## string search
  #string(SUBSTRING <string> <begin> <length> <output variable>)
  #MESSAGE( STATUS "cxx_flags: ${CMAKE_CXX_FLAGS}" )
  #string(REPLACE " " ";" FLAGS ${CMAKE_CXX_FLAGS})
  #MESSAGE( STATUS "flags: ${FLAGS}" )
  #FOREACH( flag ${FLAGS} )
  #  MESSAGE( STATUS "current flag == ${flag}" )
  #  IF( ${flag} STREQUAL "-ftemplate-depth-50" )
  #    MESSAGE( STATUS "FOUND!" )
  #  ENDIF()
  #ENDFOREACH()
  #SET( CMAKE_CXX_FLAGS "-ftemplate-depth=900 -Wall -Wno-deprecated -msse2" )
  #MESSAGE( STATUS "flags: ${CMAKE_CXX_FLAGS}")
	INCLUDE( CheckCXXCompilerFlag )
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
	IF( COMPILER_SUPPORTS_CXX11 )
		ADD_DEFINITIONS( -DCMAKE_CXX_FLAGS:STRING="${CMAKE_CXX_FLAGS} -std=c++11" )
	ELSEIF(COMPILER_SUPPORTS_CXX0X )
		ADD_DEFINITIONS( -DCMAKE_CXX_FLAGS:STRING="${CMAKE_CXX_FLAGS} -std=c++0x" )
	ELSE()
		MESSAGE(ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
	ENDIF()
  
  #MESSAGE( STATUS "flags: ${CMAKE_CXX_FLAGS}")
  
  #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ftemplate-depth=900 -fdiagnostics-color=always" )
ENDIF(UNIX) 

ExternalProject_Add( ITK
	DEPENDS ${ITK_DEPENDENCIES}
	GIT_REPOSITORY ${git_protocol}://itk.org/ITK.git #  url from where to download
	GIT_TAG v4.12.2
	SOURCE_DIR ITK
	BINARY_DIR ITK-build
	UPDATE_COMMAND ""
	PATCH_COMMAND ""
	CMAKE_GENERATOR ${gen}
	CMAKE_ARGS
		${ep_common_args}
    
		-DBUILD_EXAMPLES:BOOL=OFF # examples are not needed
		-DBUILD_SHARED_LIBS:BOOL=OFF 
		-DBUILD_TESTING:BOOL=OFF # testing the ITK build is not required
		-DITK_BUILD_ALL_MODULES:BOOL=ON
    -DITK_DYNAMIC_LOADING:BOOL=OFF
    -DModule_ITKReview:BOOL=ON
    -DVTK_DIR:PATH=${VTK_DIR}
    -DModule_ITKVtkGlue:BOOL=ON
    -DModule_LesionSizingToolkit:BOOL=ON
    -DModule_SkullStrip:BOOL=ON
    -DModule_TextureFeatures:BOOL=ON
    -DVCL_INCLUDE_CXX_0X:BOOL=ON
    -DITK_USE_SYSTEM_DCMTK:BOOL=ON
    -DModule_ITKIODCMTK:BOOL=ON
		#-DITK_LEGACY_REMOVE:BOOL=ON 
		#-DModule_ITKVideoBridgeOpenCV:BOOL=ON # [OPENCV] dependency
		#-DOpenCV_DIR:PATH=${OpenCV_DIR} # [OPENCV] dependency
		#-DModule_ITKVtkGlue:BOOL=ON # [VTK] dependency
		#-DVTK_DIR:PATH=${VTK_DIR} # [VTK] dependency
		-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} # toggle for type of build if something different that 
		-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/INSTALL
	)

SET( ITK_DIR ${CMAKE_BINARY_DIR}/ITK-build )