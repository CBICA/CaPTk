# Macro to add, link, build, and install an application and it's cwl file
MACRO( CAPTK_ADD_EXECUTABLE APPLICATION SOURCESFILES DEPENDENT_LIBS )

  #MESSAGE( AUTHOR_WARNING "APPLICATION: ${APPLICATION}" )
  #MESSAGE( AUTHOR_WARNING "SOURCESFILES: ${SOURCESFILES}" )
  #MESSAGE( AUTHOR_WARNING "DEPENDENT_LIBS: ${DEPENDENT_LIBS}" )
  
  #MESSAGE( AUTHOR_WARNING "CAPTK_BUILD_DIR: ${CAPTK_BUILD_DIR}" )
  IF( "${CAPTK_BUILD_DIR}" STREQUAL "" )
    # do nothing in this case
    #MESSAGE( AUTHOR_WARNING "CAPTK_BUILD_DIR detected as empty" )
  ELSE()
    #MESSAGE( AUTHOR_WARNING "CAPTK_BUILD_DIR detected as non-empty" )
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CAPTK_BUILD_DIR} )
    SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CAPTK_BUILD_DIR} )
    SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CAPTK_BUILD_DIR} )
  ENDIF()

  ADD_EXECUTABLE( 
    ${APPLICATION}
    ${SOURCESFILES}
  )
  
  TARGET_LINK_LIBRARIES( 
    ${APPLICATION}
    ${DEPENDENT_LIBS}
  )

  SET_TARGET_PROPERTIES( ${APPLICATION} PROPERTIES FOLDER "${CaPTk_CLIs_Folder}" )

  IF( APPLE )
    SET_TARGET_PROPERTIES( ${APPLICATION} PROPERTIES INSTALL_RPATH 
    "@loader_path/../Frameworks;@executable_path/../Frameworks;@loader_path/../lib;@executable_path/../lib;@loader_path/../../Frameworks;@executable_path/../../Frameworks;@loader_path/../../lib;@executable_path/../../lib" )
  # ELSEIF(LINUX) 
    # SET_TARGET_PROPERTIES( ${APPLICATION} PROPERTIES INSTALL_RPATH 
    # "@loader_path/../lib;@executable_path/../lib;@loader_path/../../lib;@executable_path/../../lib" )
  ENDIF()

  ADD_DEPENDENCIES( ${APPLICATION} ${LIBNAME_Applications} ${LIBNAME_FeatureExtractor} ${LIBNAME_CBICATK} )

  IF (APPLE) 
    # list (APPEND STANDALONE_APPS_LIST ${APPLICATION})
    INSTALL( 
      TARGETS ${APPLICATION}
      BUNDLE DESTINATION .
      RUNTIME DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
      LIBRARY DESTINATION ${EXE_NAME}.app/Contents/Resources/lib
      CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    ) 

  ELSE()
    INSTALL( 
      TARGETS ${APPLICATION}
      BUNDLE DESTINATION .
      RUNTIME DESTINATION bin
      LIBRARY DESTINATION lib
      CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

  ENDIF()

  # Add test for run tests
  #ADD_TEST( NAME ${APPLICATION}_rt COMMAND ${APPLICATION} -rt )

  IF (NOT CAPTK_DISABLE_CWL)
    CWL_INSTALL(${APPLICATION})
  ENDIF()
  
ENDMACRO()

# Macro to generate and install a cwl file after target application is built
MACRO(CWL_INSTALL APPLICATION)

  # Post build cwl generation
  add_custom_command(
    TARGET ${APPLICATION}
    POST_BUILD
    COMMAND ${APPLICATION} -cwl
    COMMENT "Generating cwl file for ${APPLICATION}..."
  )
  
  # Post build cwl generation
  add_custom_command(
    TARGET ${APPLICATION}
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CAPTK_BUILD_DIR}/${CMAKE_CONFIGURATION_TYPES}/${APPLICATION}.cwl ${DATA_DIR}/cwlFiles/${APPLICATION}.cwl
    COMMENT "Updating cwl File for ${APPLICATION}..."
  )
  
  #IF (APPLE) 
  #  list (APPEND STANDALONE_APPS_LIST ${APPLICATION})
  #  INSTALL( 
  #    FILES ${PROJECT_BINARY_DIR}/${APPLICATION}.cwl
  #    DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
  #  ) 
  #
  #ELSEIF(WIN32)
  #  INSTALL( 
  #    FILES ${PROJECT_BINARY_DIR}/$<CONFIGURATION>/${APPLICATION}.cwl
  #    DESTINATION bin
  #  )
  #
  #ELSE()
  #  INSTALL( 
  #    FILES ${PROJECT_BINARY_DIR}/${APPLICATION}.cwl
  #    DESTINATION bin
  #  )
  #
  #ENDIF()

ENDMACRO()

# macro to find all sub-directories
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children
    RELATIVE ${curdir} ${curdir}/*
    PATTERN "svn" EXCLUDE
  )
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

# a common version number is always maintained in CaPTk and all its applications
MACRO( CAPTK_ADD_PROJECT NEW_PROJECT_NAME NEW_PROJECT_VERSION )
  
  IF( "${PROJECT_VERSION}" STREQUAL "" )
    # this basically means that packaging will not be happening since the project is independent of CaPTk 
    PROJECT( ${NEW_PROJECT_NAME} )
    SET( PROJECT_VERSION "${NEW_PROJECT_VERSION}" )
    ADD_DEFINITIONS( -DPROJECT_VERSION="${PROJECT_VERSION}" )
    SET( CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Default configuration types" FORCE )
  ELSE()
    PROJECT( ${NEW_PROJECT_NAME} )  
  ENDIF()  
  
ENDMACRO()

# macro to handle initial setup of projects (no library dependency management)
MACRO( CAPTK_INITIAL_SETUP )

CMAKE_POLICY(SET CMP0048 NEW) # only with CMake 3
CMAKE_POLICY(SET CMP0025 NEW) # only with CMake 3

FIND_PACKAGE( ITK )
IF( EXISTS ${ITK_USE_FILE} )
  INCLUDE( ${ITK_USE_FILE} )
ENDIF()

SET(CMAKE_CXX_STANDARD 11 CACHE STRING "C++ standard to be used" )
SET(CMAKE_CXX_STANDARD_REQUIRED YES) 
SET_PROPERTY( GLOBAL PROPERTY USE_FOLDERS ON )

SET( CACHED_INCLUDE_DIRS
  ${CACHED_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/../common_includes/
  ${PROJECT_SOURCE_DIR}/../../common/inc
  ${PROJECT_SOURCE_DIR}/src/
  ${PROJECT_SOURCE_DIR}/src/depends/
  CACHE STRING "All include directories" FORCE
)
#MESSAGE( STATUS "[DEBUG] CACHED_INCLUDE_DIRS@Macro: ${CACHED_INCLUDE_DIRS}" )

FILE( GLOB_RECURSE CURRENT_APPLICATION_DEPENDS "${PROJECT_SOURCE_DIR}/src/depends/*" )
FILE( GLOB CURRENT_APPLICATION_SOURCES "${PROJECT_SOURCE_DIR}/src/${PROJECT_NAME}*" )

SET( APPLICATION_SOURCES ${CURRENT_APPLICATION_SOURCES} CACHE STRING "All application sources" FORCE )
SET( APPLICATION_DEPENDS ${CURRENT_APPLICATION_DEPENDS} CACHE STRING "All application depends" FORCE )

# set the install path
IF( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  SET( CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/install CACHE PATH "Changed the default install path" FORCE  )
ENDIF( )

#  set cmake up
IF(WIN32)
  SET(CMAKE_C_USE_RESPONSE_FILE_FOR_OBJECTS 1)
  SET(CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS 1)
  SET(CMAKE_C_USE_RESPONSE_FILE_FOR_INCLUDES 1)
  SET(CMAKE_CXX_USE_RESPONSE_FILE_FOR_INCLUDES 1)
  SET(CMAKE_C_RESPONSE_FILE_LINK_FLAG "@")
  SET(CMAKE_CXX_RESPONSE_FILE_LINK_FLAG "@")
  SET(CMAKE_C_CREATE_STATIC_LIBRARY "<CMAKE_AR> rc <TARGET> <LINK_FLAGS> <OBJECTS>") 
  SET(CMAKE_CXX_CREATE_STATIC_LIBRARY "<CMAKE_AR> rc <TARGET> <LINK_FLAGS> <OBJECTS>") 
ENDIF()

IF(NOT WIN32)
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
ENDIF()

# find openmp flags
IF(APPLE)
  SET(OPENMP_LIBRARIES "${CMAKE_C_COMPILER}/../../lib")
  SET(OPENMP_INCLUDES "${CMAKE_C_COMPILER}/../../include")
  
  # MESSAGE( AUTHOR_WARNING "Apple_C_Compiler: ${CMAKE_C_COMPILER}")
  
  SET(OpenMP_C "${CMAKE_C_COMPILER}")
  SET(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
  SET(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
  SET(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
  SET(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
  SET(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
  SET(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
  SET(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
  SET(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
  SET(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  SET(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  SET(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  
  INCLUDE_DIRECTORIES("${OPENMP_INCLUDES}")
  LINK_DIRECTORIES("${OPENMP_LIBRARIES}")

 ENDIF()

FIND_PACKAGE(OpenMP REQUIRED)
SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

#MESSAGE(AUTHOR_WARNING "OpenMP is enabled, you may see issues on Azure")

# add a few compiler options for GCC and MSVC
IF( CMAKE_COMPILER_IS_GNUCXX )
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -O3" )
ENDIF()

if(APPLE)
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-dev -Wno-narrowing -Wall" )
endif()
  
ENDMACRO()