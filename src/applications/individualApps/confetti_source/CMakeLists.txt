CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#PROJECT(Confetti)
#SET( PROJECT_VERSION_MAJOR 1)
#SET( PROJECT_VERSION_MINOR 0)
#SET( PROJECT_VERSION_PATCH 0)
#SET( PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")


FUNCTION(FIND_PYTHON_MODULE module)
  #reference http://www.cmake.org/pipermail/cmake/2011-January/041666.html
  include(FindPackageHandleStandardArgs)
  string(TOUPPER ${module} module_upper)
  if(NOT PY_${module_upper})
    if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
      set(${module}_FIND_REQUIRED TRUE)
    endif()
    # A module's location is usually a directory, but for binary modules
    # it's a .so file.
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" 
      "import re, ${module}; print(re.compile('/__init__.py.*').sub('',${module}.__file__))"
      RESULT_VARIABLE _${module}_status 
      OUTPUT_VARIABLE _${module}_location
      ERROR_QUIET 
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT _${module}_status)
      set(PY_${module_upper} ${_${module}_location} CACHE STRING 
        "Location of Python module ${module}")
    endif(NOT _${module}_status)
  endif(NOT PY_${module_upper})
  find_package_handle_standard_args(PY_${module} DEFAULT_MSG PY_${module_upper})
ENDFUNCTION(FIND_PYTHON_MODULE)

FIND_PACKAGE( PythonInterp 2.7 REQUIRED )
SET(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} CACHE FILEPATH "python interpreter")
FIND_PYTHON_MODULE(numpy REQUIRED)
FIND_PYTHON_MODULE(vtk REQUIRED)
FIND_PYTHON_MODULE(PyQt4 REQUIRED)
ADD_DEFINITIONS(-DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE})
IF(PY_NUMPY AND PY_VTK AND PY_PYQT4)
  MESSAGE("found numpy, vtk and pyQT4 for ConfettiGUI")
  #ADD_SUBDIRECTORY(${SUBDIRPATH})
ELSE()
  MESSAGE(WARNING "Canot find python packages required for ConfettiGUI(numpy, vtk and pyQT4). Please ensure Python dependencies are met (see ${captk_source}/src/applications/binaries/confetti/pyGUI/requirements for details).")
ENDIF()
      

IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-function")
ENDIF() 
CMAKE_POLICY(SET CMP0003 NEW)
FIND_PACKAGE(OpenCV REQUIRED)
FIND_PACKAGE(ITK REQUIRED)
INCLUDE("${ITK_USE_FILE}")
SET(ConfettiSrcs src/main.cpp
	src/confetti.cpp)
SET(ConfettiHeaders 
  src/confetti.h
  src/dtiUtils.h
  src/multinomialModel.h)

ADD_EXECUTABLE(Confetti ${ConfettiSrcs} ${ConfettiHeaders})
TARGET_LINK_LIBRARIES(Confetti ${OpenCV_LIBS} ${ITK_LIBRARIES})

add_library(ConfettiCore SHARED
	src/confetti.cpp
	src/confetti.h
	src/multinomialModel.h
	src/dtiUtils.h
	src/confettiApi.cpp
	)
TARGET_LINK_LIBRARIES(ConfettiCore ${OpenCV_LIBS} ${ITK_LIBRARIES})

# install all targets and scripts
INSTALL( TARGETS Confetti 
  DESTINATION bin  
  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

INSTALL( TARGETS ConfettiCore 
  DESTINATION lib  
  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

INSTALL( FILES pyGUI/ConfettiCore.py
  DESTINATION libexec  
  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

INSTALL( FILES pyGUI/ConfettiGUI.py
  DESTINATION bin  
  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

SET( DATA_DIR "${PROJECT_SOURCE_DIR}/data" )
#INCLUDE(CTest)
#ENABLE_TESTING()
#ADD_TEST( NAME Connectivity COMMAND Confetti -test Connectivity)
#ADD_TEST( NAME Clustering COMMAND Confetti -test Clustering)
#ADD_TEST( NAME Extraction COMMAND Confetti -test Extraction)
