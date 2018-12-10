# This CMakeLists file will add the folders "GUI", "CBICA_Toolkit" and "Applications" to the compilation

# It will also assess the number of Applications present (using "file(GLOB_RECURSE HDRS *.h)")

# Add a variable to the local cache which denotes the number of applications, so that it can be accessed by the c++ source code and add the requisite number of menu items

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/applications)

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/cbica_toolkit)

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/gui)

ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/common)

IF(WIN32)
	#ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/lib/Mathworks/MCR_R2013A_Win64)
	#ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/src/lib/ITK-SNAP/3.4_Win64)
ENDIF(WIN32)

SET( EXE_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/CAPTk.cpp"
  CACHE STRING "Executable sources"
)
