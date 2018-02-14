## The Superbuild script is used to automatically download and build project dependencies. 

# Using GIT to download third party libraries. An SVN/BitBucket URL will also work the same way
FIND_PACKAGE( Git REQUIRED )

OPTION( USE_GIT_PROTOCOL "If behind a firewall turn this off to use https instead." ON)

SET(git_protocol "git")
IF(NOT USE_GIT_PROTOCOL)
	SET(git_protocol "https")
ENDIF()

SET( CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
)

INCLUDE( ExternalProject )

# Compute -G arg for configuring external projects with the same CMake generator:
IF(CMAKE_EXTRA_GENERATOR)
	SET(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
ELSE()
	SET(gen "${CMAKE_GENERATOR}" )
ENDIF()

SET( NewCore_DEPENDENCIES )

# Automatic VTK build and link
#OPTION( USE_VTK "Build VTK v5.10.1" OFF )
#IF( ${USE_VTK} )
#	FIND_PACKAGE( VTK REQUIRED )
#	IF( NOT VTK_DIR )
#		MESSAGE( STATUS "VTK not found on system. Building from source..." )
#		INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/External-VTK.cmake )
#	ENDIF(  )
#	SET( NewCore_DEPENDENCIES ${NewCore_DEPENDENCIES} VTK )
#ENDIF()

# Automatic DCMTK build and link
OPTION( USE_DCMTK "Build DCMTK v3.6.1" OFF )
IF( ${USE_DCMTK} )
	FIND_PACKAGE( DCMTK REQUIRED )
	IF( NOT DCMTK_DIR )
		MESSAGE( STATUS "DCMTK not found on system. Building from source..." )
		INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/External-DCMTK.cmake )
	ENDIF(  )
	SET( NewCore_DEPENDENCIES ${NewCore_DEPENDENCIES} DCMTK )
ENDIF()

# Automatic ITK build and link
OPTION( USE_ITK "Use ITK v4.12.2" ON )
IF( ${USE_ITK} )
	FIND_PACKAGE( ITK REQUIRED )
	IF( NOT ITK_DIR )
		MESSAGE( STATUS "ITK not found on system. Building from source..." )
		INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/External-ITK.cmake )
	ENDIF( )
	SET( NewCore_DEPENDENCIES ${NewCore_DEPENDENCIES} ITK )
ENDIF( )
