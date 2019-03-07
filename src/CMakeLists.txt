# This CMakeLists file will add the folders "GUI", "CBICA_Toolkit" and "Applications" to the compilation

# It will also assess the number of Applications present (using "file(GLOB_RECURSE HDRS *.h)")

# Add a variable to the local cache which denotes the number of applications, so that it can be accessed by the c++ source code and add the requisite number of menu items

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/cbica_toolkit)

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/applications)

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/common)

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/view)

IF(WIN32)
	#ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/lib/Mathworks/MCR_R2013A_Win64)
	#ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/lib/ITK-SNAP/3.4_Win64)
ENDIF(WIN32)

SET( EXE_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/CaPTk.cxx"
  CACHE STRING "Executable sources"
)
