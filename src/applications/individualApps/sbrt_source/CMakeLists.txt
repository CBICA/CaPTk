#CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#
#PROJECT( SBRT_Lung )

SET( SBRT_PROJECT_NAME SBRT_Lung )
SET ( Executable_1 ${SBRT_PROJECT_NAME}_Segment )
SET ( Executable_2 ${SBRT_PROJECT_NAME}_Analyze )

#SET(PROJECT_VERSION "1.0.0")
#ADD_DEFINITIONS(-DPROJECT_VERSION="${PROJECT_VERSION}")
 
#Find libraries
FIND_PACKAGE( ITK REQUIRED )
INCLUDE( ${ITK_USE_FILE} )
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
FIND_PACKAGE(OpenCV  REQUIRED )
include_directories(${OpenCV_INCLUDE_DIRS})

IF (NOT APPLE) {
  FIND_PACKAGE(OpenMP REQUIRED)
  SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
}

#INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/src/segment/ )
#INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/src/analyze/ )
#INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/)

SET( segment_src 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/segment/Joint_Segm.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/segment/Joint_Segm.cpp
)

SET( analyze_src 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/analyze/LBPFeatures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/analyze/LBPFeatures.cpp  
  ${CMAKE_CURRENT_SOURCE_DIR}/src/analyze/generateTextureFeatures.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/analyze/generateTextureFeatures.cpp

)


SET( Common_dependents

  #${CMAKE_CURRENT_SOURCE_DIR}/common_includes/cbicaPreProcessImage.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaLogging.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaLogging.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaUtilities.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaUtilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaCmdParser.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaCmdParser.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaITKImageInfo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaITKImageInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common_includes/cbicaITKSafeImageIO.h

)
### Building the documentation
## The following 4 lines are for generation of public documentation
#SET( BUILD_DOCUMENTATION ON "Do you want to build Documentation (needs Doxygen)?" )
#MARK_AS_ADVANCED(BUILD_DOCUMENTATION)
#IF( BUILD_DOCUMENTATION )
# ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/docs)
#ENDIF ()

# add c++11 flag if it is GCC v<5
#IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#  IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
#    INCLUDE( CheckCXXCompilerFlag )
#    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#    IF( COMPILER_SUPPORTS_CXX11 )
#	MESSAGE(STATUS "Compiler ${CMAKE_CXX_COMPILER} version ${CMAKE_CXX_COMPILER_VERSION} supports C++11")
#    	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#    elseif(COMPILER_SUPPORTS_CXX0X )
#	MESSAGE(STATUS "Compiler ${CMAKE_CXX_COMPILER} version ${CMAKE_CXX_COMPILER_VERSION} supports C++0x")
#    	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#    else()
#    	MESSAGE(ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#    endif()
#  ELSE(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
#	MESSAGE(STATUS "Compiler ${CMAKE_CXX_COMPILER} version ${CMAKE_CXX_COMPILER_VERSION} is above 5.0, known support for C++11, C++0x")
#    	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#    	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#  ENDIF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
#ENDIF() 


#set(BUILD_TESTING TRUE CACHE BOOL "Build testing")
#if(BUILD_TESTING)
#	enable_testing()
#	ENABLE_TESTING()
#	add_subdirectory(testing)
#	INCLUDE(CTest)
#endif(BUILD_TESTING)


SET( CMD_MODE OFF CACHE BOOL "Builds exe suitable for put in to CaPTk." )
ADD_DEFINITIONS(-DCMD_MODE)



ADD_EXECUTABLE(
  ${Executable_1} 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/segment/main.cxx
  ${segment_src}
  ${Common_dependents}
)
ADD_EXECUTABLE(
  ${Executable_2} 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/analyze/main.cxx
  ${analyze_src}
  ${Common_dependents}
)

# Link the libraries to be used
IF( ITKVtkGlue_LOADED )
  TARGET_LINK_LIBRARIES(
    ${Executable_1}
  	${Glue}  
    ${VTK_LIBRARIES} 
    ${ITK_LIBRARIES}
    ${OpenCV_LIBS}
  )
  TARGET_LINK_LIBRARIES(
    ${Executable_2}
  	${Glue}  
    ${VTK_LIBRARIES} 
    ${ITK_LIBRARIES}
    ${OpenCV_LIBS}
  )
ELSE()
  TARGET_LINK_LIBRARIES(
    ${Executable_1}
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES} 
    ${OpenCV_LIBS}
  )
  TARGET_LINK_LIBRARIES(
    ${Executable_2}
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES} 
    ${OpenCV_LIBS}
  )
ENDIF()

INSTALL( TARGETS ${Executable_1} 
  DESTINATION bin  
  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
INSTALL( TARGETS ${Executable_2} 
  DESTINATION bin  
  CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)