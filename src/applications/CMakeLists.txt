# This CMakeLists file will add the source files in this directory for compilation

# common_includes will have all the functions and classes from CBICA_Toolkit_ITK branch and more

# each header file, "$application.h" : main entry point of an application whose input will be either a single "itk::Image" or a "std::vector< itk::Image >"

# each application-specific include folder, "/$application_includes" : dependencies which are not taken care by the default dependencies (i.e., ITK and OpenCV), for example libSVM for recurrence 

FILE( GLOB 
  COMMON_DEPENDENCIES
  "${CMAKE_CURRENT_SOURCE_DIR}/common_includes/*.*"

)
SET( SOURCES "" )
SET( DEPENDENCIES ${COMMON_DEPENDENCIES} )
SET( INCLUDES "")

SET( APPLICATIONS
  WhiteStripe
  GeodesicSegmentation
  EGFRvIIISurrogateIndex
  RecurrenceEstimator
  SurvivalPredictor
  MolecularSubtypePredictor
  ImagingSubtypePredictor
  PopulationAtlases
  DirectionalityEstimate
  FeatureExtraction
)
SET( CAPTK_APP_LIST "" )

SET( CAPTK_BRAIN_APPS "" )
SET( CAPTK_BREAST_APPS "" )
SET( CAPTK_LUNG_APPS "" )
SET( CAPTK_MISC_APPS "" )

SET(PREPROCESSING_ALGOS "DCM2NIfTI Denoise-SUSAN Register BiasCorrect-N3 HistogramMatching SkullStripping PCA PerfusionDerivatives DiffusionDerivatives" CACHE STRING "All Preprocessing algorithms available")

## pass the application under consideration and the specific APP_LIST (either BRAIN, BREAST, LUNG or MISC that you need to add it to, if found)
#MACRO( ADD_APPLICATION application APP_LIST )
#  # ensure that specified application has been found
#  IF( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${application}.h")
#  
#    SET( SOURCES ${SOURCES} "${CMAKE_CURRENT_SOURCE_DIR}/${application}.h" )
#  
#  # if implementation is present, include that for compilation
#    IF( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${application}.cpp")
#      SET( SOURCES ${SOURCES} "${CMAKE_CURRENT_SOURCE_DIR}/${application}.cpp" )
#    ENDIF()
#  
#    IF( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/" )
#      # include dependency directory of ${application} for compilation
#      SET( INCLUDES ${INCLUDES} "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/" )
#    
#      # add all dependencies of ${application} to compilation
#      FILE( GLOB 
#        APPLICATION_SPECIFIC_DEPENDENCIES
#        "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/*.h"
#        "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/*.cpp"
#      )
#    ENDIF()
#  
#    # keep storing all dependencies
#    SET( DEPENDENCIES ${DEPENDENCIES} ${APPLICATION_SPECIFIC_DEPENDENCIES} )
#  
#    # add to variable which will be used in source code to generate menu items
#    SET( ${APP_LIST} "${${APP_LIST}} ${application}" )
#    
#  ENDIF()
#ENDMACRO()

MACRO( ADD_APPLICATION application)
  # ensure that specified application has been found
  IF( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${application}.h")
  
    SET( SOURCES ${SOURCES} "${CMAKE_CURRENT_SOURCE_DIR}/${application}.h" )
  
  # if implementation is present, include that for compilation
    IF( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${application}.cpp")
      SET( SOURCES ${SOURCES} "${CMAKE_CURRENT_SOURCE_DIR}/${application}.cpp" )
    ENDIF()
  
    IF( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/" )
      # include dependency directory of ${application} for compilation
      SET( INCLUDES ${INCLUDES} "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/" )
    
      # add all dependencies of ${application} to compilation
      FILE( GLOB 
        APPLICATION_SPECIFIC_DEPENDENCIES
        "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/${application}_includes/*.cpp"
      )
    ENDIF()
  
    # keep storing all dependencies
    SET( DEPENDENCIES ${DEPENDENCIES} ${APPLICATION_SPECIFIC_DEPENDENCIES} )
  
    # add to variable which will be used in source code to generate menu items
    SET( CAPTK_APP_LIST "${CAPTK_APP_LIST} ${application}" )
    
  ENDIF()
ENDMACRO()

FOREACH(application ${APPLICATIONS})
  
  IF ( BUILD_EGFRvIII AND (${application} STREQUAL EGFRvIIISurrogateIndex))
    ADD_APPLICATION( ${application} )
  ENDIF()
  
  IF ( BUILD_RECURRENCE AND (${application} STREQUAL RecurrenceEstimator))
    ADD_APPLICATION( ${application} )
  ENDIF()
  
  IF ( BUILD_GEODESIC AND (${application} STREQUAL GeodesicSegmentation))
    ADD_APPLICATION( ${application} )
  ENDIF()
  
  IF ( BUILD_SURVIVAL AND (${application} STREQUAL SurvivalPredictor))
    ADD_APPLICATION( ${application} )
  ENDIF()
  
  IF ( BUILD_ISUBTYPE AND (${application} STREQUAL ImagingSubtypePredictor))
    ADD_APPLICATION( ${application} )
  ENDIF()
  
 
  IF ( BUILD_MSUBTYPE AND (${application} STREQUAL MolecularSubtypePredictor))
    ADD_APPLICATION( ${application} )
  ENDIF()

    IF ( BUILD_ATLAS AND (${application} STREQUAL PopulationAtlases))
    ADD_APPLICATION( ${application} )
  ENDIF()

  IF ( (${application} STREQUAL DirectionalityEstimate))
    ADD_APPLICATION( ${application} )
  ENDIF()
  
  IF( (${application} STREQUAL WhiteStripe) )
    ADD_APPLICATION( ${application} )
  ENDIF()
  
  IF ( (${application} STREQUAL FeatureExtraction))
    SET( CAPTK_APP_LIST "${CAPTK_APP_LIST} ${application}" )
  ENDIF()
  
ENDFOREACH()

## this is for all stand-alone applications which do not integrate at source code level with CaPTk
SET( CAPTK_APP_LIST_PY_GUI "" )

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

## download and extract the prepared binaries of external apps

    ## download exe from url
    #IF (PACKAGE_PROJECT)
    #  FILE(DOWNLOAD "${URL_GOES_HERE}" "${SUBDIRPATH}/libra.zip" TIMEOUT 600 STATUS STATUS_CODE)
    #  IF(NOT STATUS_CODE EQUAL 0)
    #    message(FATAL_ERROR "Failed to download LIBRA. Status=${STATUS_CODE}")
    #  ENDIF()
    #ENDIF()
    #EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xz "${SUBDIRPATH}/libra.zip"
    #              WORKING_DIRECTORY "${SUBDIRPATH}/LIBRA"
    #              RESULT_VARIABLE RESULT_CODE)
    #IF(NOT RESULT_CODE EQUAL 0)
    #  message(FATAL_ERROR "Failed to unpack LIBRA from ${SUBDIRPATH}/libra.zip to ${SUBDIRPATH}/LIBRA with error ${RESULT_CODE}")
    #ENDIF()

SET( FILENAME_TO_EXTRACT "" )
SET( DOWNLOAD_URL "" )
SET( DOWNLOADED_APPS_DIR ${PROJECT_BINARY_DIR}/externalApps )
SET( SOURCE_APPLICATIONS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/individualApps )

# putting this condition so that the extraction doesn't happen every single time
IF( NOT EXISTS "${DOWNLOADED_APPS_DIR}" )  

  IF(WIN32)
    SET( FILENAME_TO_EXTRACT "binaries_windows")
    SET( DOWNLOAD_URL "ftp://www.nitrc.org/var/home/groups/captk/downloads/sourceBuild_extra_windows.zip")
  ELSE()
    SET( FILENAME_TO_EXTRACT "binaries_linux")
    SET( DOWNLOAD_URL "ftp://www.nitrc.org/var/home/groups/captk/downloads/sourceBuild_extra_linux.zip")
  ENDIF()

  SET( FILE_TO_EXTRACT ${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip)

  # download exe from url
  MESSAGE( STATUS "Downloading pre-compiled external applications" )
  FILE(DOWNLOAD "${DOWNLOAD_URL}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
  IF(NOT STATUS_CODE EQUAL 0)
    message(FATAL_ERROR "Failed to download Precompiled packages. Status=${STATUS_CODE}")
  ENDIF()

  FILE(MAKE_DIRECTORY ${DOWNLOADED_APPS_DIR})

  MESSAGE( STATUS "Extracting pre-compiled external applications" )
  IF( EXISTS "${FILE_TO_EXTRACT}" )

    EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
      WORKING_DIRECTORY ${DOWNLOADED_APPS_DIR}
      RESULT_VARIABLE RESULT_CODE
    )
    
    IF(NOT RESULT_CODE EQUAL 0)
      MESSAGE( WARNING "Extracting the pre-compiled applications failed" )
    ENDIF()
    
    SUBDIRLIST(EXTERNAL_APPS ${DOWNLOADED_APPS_DIR}/${FILENAME_TO_EXTRACT})
    
    # copy the prepared app into the appropriate folder
    FOREACH(external_app ${EXTERNAL_APPS})
      
      FILE(MAKE_DIRECTORY ${SOURCE_APPLICATIONS_PATH}/${external_app})
      
      EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${DOWNLOADED_APPS_DIR}/${FILENAME_TO_EXTRACT}/${external_app} ${SOURCE_APPLICATIONS_PATH}/${external_app}
      )
    
    ENDFOREACH()
    
  ENDIF()
ENDIF()

SUBDIRLIST(SUBDIRECTORIES ${SOURCE_APPLICATIONS_PATH})
INCLUDE(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

FOREACH(subdir ${SUBDIRECTORIES})
   
  SET(SUBDIRPATH "${SOURCE_APPLICATIONS_PATH}/${subdir}")
  
  IF (BUILD_LIBRASINGLE AND (${subdir} STREQUAL "libra"))
    #FIND_PACKAGE(MATLAB REQUIRED)
    #MATLAB_GET_RELEASE_NAME_FROM_VERSION(MATLAB_VERSION)
    #MESSAGE(STATUS "[debug]matlab_Version=${MATLAB_VERSION}")
    
    #TBD: discuss name with Mark
    SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} librasingle" CACHE STRING "Available stand-alone apps" FORCE)
    # ALSO add this to CAPTK_BRREAST_APPS    
    
    IF(WIN32)
      INSTALL(FILES ${SUBDIRPATH}/libra.exe
        DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
      INSTALL(FILES ${SUBDIRPATH}/libra.bat
        DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )  
      INSTALL(DIRECTORY "${SUBDIRPATH}/Mathworks" DESTINATION bin)
    ELSE()
      INSTALL(FILES ${SUBDIRPATH}/bin/libra
        DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
      INSTALL(FILES ${SUBDIRPATH}/libexec/libra
        DESTINATION libexec
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
      INSTALL(FILES ${SUBDIRPATH}/MCRInstaller.zip
        DESTINATION preinstall
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
    ENDIF()
  
  ENDIF()  
    
  IF(BUILD_LIBRABATCH AND (${subdir} STREQUAL "libra"))
    #FIND_PACKAGE(MATLAB REQUIRED)
    #MATLAB_GET_RELEASE_NAME_FROM_VERSION(MATLAB_VERSION)
    #MESSAGE(STATUS "[debug]matlab_Version=${MATLAB_VERSION}")
    
    ## TBD: this needs to be discussed with Mark
    SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} librabatch" CACHE STRING "Available stand-alone apps" FORCE)
    # ALSO add this to CAPTK_BRREAST_APPS
    
    IF(WIN32)      
      INSTALL(FILES ${SUBDIRPATH}/libra.exe
        DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
      INSTALL(FILES ${SUBDIRPATH}/libra.bat
        DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )  
      INSTALL(DIRECTORY "${SUBDIRPATH}/Mathworks" DESTINATION bin)
    ELSE()    
      INSTALL(FILES ${SUBDIRPATH}/bin/libra
        DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
      INSTALL(FILES ${SUBDIRPATH}/libexec/libra
        DESTINATION libexec
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
      INSTALL(FILES ${SUBDIRPATH}/MCRInstaller.zip
        DESTINATION preinstall
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
    ENDIF()
    
  ENDIF()
  
  IF (BUILD_ITKSNAP AND (${subdir} STREQUAL "itksnap"))
  
    SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} itksnap" CACHE STRING "Available stand-alone apps" FORCE )
    # ALSO add this to CAPTK_BRAIN_APPS
  
    IF (WIN32)
      FILE(GLOB allFiles ${SUBDIRPATH}/bin/*.*)
      INSTALL(FILES ${allFiles}
        DESTINATION bin
        #FILES_MATCHING
        #PATTERN "*.*"
        #PATTERN ".svn" EXCLUDE
      )
      INSTALL(DIRECTORY ${SUBDIRPATH}/bin/platforms
        DESTINATION bin
        #FILES_MATCHING
        #PATTERN "*.*"
        #PATTERN ".svn" EXCLUDE
      )
    ELSE()  
      INSTALL( FILES ${SUBDIRPATH}/bin/itksnap
        DESTINATION bin
        #FILES_MATCHING
        #PATTERN "*.*"
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE

      )
      #FILE(GLOB allFiles ${SUBDIRPATH}/lib/snap-3.6.0-rc1/*.*)
      #INSTALL(FILES ${allFiles}
      #  DESTINATION lib/snap-3.6.0-rc1
        #FILES_MATCHING
        #PATTERN "*.*"
        #PATTERN ".svn" EXCLUDE
      #  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      #)
      INSTALL(DIRECTORY ${SUBDIRPATH}/lib/snap-3.6.0-rc1
        DESTINATION lib
        USE_SOURCE_PERMISSIONS
        #FILES_MATCHING
        #PATTERN "*.*"
        #PATTERN ".svn" EXCLUDE
        #PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
      )
    ENDIF()
    
    ## download exe from url
    #IF (PACKAGE_PROJECT)
    #  FILE(DOWNLOAD "${URL_GOES_HERE}" "${SUBDIRPATH}/ITK-SNAP.zip" TIMEOUT 600 STATUS STATUS_CODE)
    #  IF(NOT STATUS_CODE EQUAL 0)
    #    message(FATAL_ERROR "Failed to download ITK-SNAP. Status=${STATUS_CODE}")
    #  ENDIF()
    #ENDIF()
    #EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xz "${SUBDIRPATH}/ITK-SNAP.zip"
    #              WORKING_DIRECTORY "${SUBDIRPATH}/ITK-SNAP"
    #              RESULT_VARIABLE RESULT_CODE)
    #IF(NOT RESULT_CODE EQUAL 0)
    #  message(FATAL_ERROR "Failed to unpack ITK-SNAP from ${SUBDIRPATH}/ITK-SNAP.zip to ${SUBDIRPATH}/ITK-SNAP with error ${RESULT_CODE}")
    #ENDIF()
    
  ENDIF()
  
  #IF (BUILD_WHITESTRIPE AND (${subdir} STREQUAL "whitestripe_source")) # this is not needed since WS has been pulled in as a source-integrated app
  #  
  #  #ADD_SUBDIRECTORY(${SUBDIRPATH})
  #  #MESSAGE( "External project - WhiteStripe" )
  #  #ExternalProject_Add( ${subdir}
  #  #  DOWNLOAD_COMMAND ""
  #  #  PATCH_COMMAND ""
  #  #  SOURCE_DIR ${SUBDIRPATH}
  #  #  #DEPENDS ${ITK_DIR}; ${OpenCV_DIR}
  #  #  #PREFIX ${PROJECT_BINARY_DIR}/${subdir}
  #  #  BINARY_DIR ${PROJECT_BINARY_DIR}/${subdir}
  #  #  INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  #  #  UPDATE_COMMAND ""
  #  #  PATCH_COMMAND ""
  #  #  CMAKE_GENERATOR ${gen}
  #  #  CMAKE_ARGS
  #  #    ${ep_common_args}
  #  #    -DITK_DIR=${ITK_DIR}
  #  #    -DOpenCV_DIR=${OpenCV_DIR}
  #  #    -DBUILD_TESTING:BOOL=OFF # testing the ITK build is not required
  #  #    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} # toggle for type of build if something different that 
  #  #)
  #
  #  SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} WhiteStripe" CACHE STRING "Available stand-alone apps" FORCE)
  #
  #ENDIF()
  
  IF (BUILD_CONFETTI AND (${subdir} STREQUAL "confetti"))  
  
    SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} confetti" CACHE STRING "Available stand-alone apps" FORCE)
    # ALSO add this to CAPTK_BRAIN_APPS
    IF (${subdir} STREQUAL "confetti") # this is only valid for WIN32 since Confetti is distributed as a compiled EXE for it
      
      FILE(GLOB allFiles ${SUBDIRPATH}/*.*)
      INSTALL(FILES ${allFiles}
        DESTINATION bin
        #FILES_MATCHING
        #PATTERN "*.*"
        #PATTERN ".svn" EXCLUDE
      )
      INSTALL(DIRECTORY ${SUBDIRPATH}/ConfettiBundleAtlas
        DESTINATION bin
        #FILES_MATCHING
        #PATTERN "*.*"
        #PATTERN ".svn" EXCLUDE
      )
      
    ENDIF()
    
  ENDIF()
    
  IF (BUILD_SBRT AND (${subdir} STREQUAL "sbrt_source"))
    
    #ADD_SUBDIRECTORY(${SUBDIRPATH})
    #MESSAGE( "External project - SBRT" )
    #ExternalProject_Add( ${subdir}
    #  DOWNLOAD_COMMAND ""
    #  PATCH_COMMAND ""
    #  SOURCE_DIR ${SUBDIRPATH}
    #  #DEPENDS ${ITK_DIR}${OpenCV_DIR}
    #  #PREFIX ${PROJECT_BINARY_DIR}/${subdir}
    #  BINARY_DIR ${PROJECT_BINARY_DIR}/${subdir}
    #  UPDATE_COMMAND ""
    #  PATCH_COMMAND ""
    #  CMAKE_GENERATOR ${gen}
    #  INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    #  CMAKE_ARGS
    #    ${ep_common_args}
    #    -DITK_DIR=${ITK_DIR}
    #    -DOpenCV_DIR=${OpenCV_DIR}
    #    -DBUILD_TESTING:BOOL=OFF # testing the ITK build is not required
    #    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} # toggle for type of build if something different that 
    #)

    SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} SBRT_Segment SBRT_Analyze" CACHE STRING "Available stand-alone apps" FORCE)
    # ALSO add this to CAPTK_LUNG_APPS
  
  ENDIF()
  
  IF (${subdir} STREQUAL "dcm2nii") # this is a pre-processing algorithm and should be present as a pre-compiled binary package
  
    SET( DCM2NII_EXE "" )
    IF (WIN32)
      SET( DCM2NII_EXE "${SUBDIRPATH}/dcm2nii.exe" )
    ELSEIF (UNIX)
      SET( DCM2NII_EXE "${SUBDIRPATH}/dcm2nii" )
    ENDIF()
    
    INSTALL(FILES ${DCM2NII_EXE}
       DESTINATION bin
       PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
     )
     
  ENDIF()
  
ENDFOREACH()
  

SET( CAPTK_APP_LIST ${CAPTK_APP_LIST} CACHE STRING "Available Code-level Applications" FORCE)
#SET( CAPTK_BRAIN_APPS ${CAPTK_BRAIN_APPS} CACHE STRING "Available Brain Applications" FORCE )
#SET( CAPTK_BREAST_APPS ${CAPTK_BREAST_APPS} CACHE STRING "Available Breast Applications" FORCE )
#SET( CAPTK_LUNG_APPS ${CAPTK_LUNG_APPS} CACHE STRING "Available Lung Applications" FORCE )
#SET( CAPTK_MISC_APPS ${CAPTK_MISC_APPS} CACHE STRING "Available Misc. Applications" FORCE )
SET( APPLICATION_SOURCES  "${SOURCES}" CACHE STRING "Code-level Application sources" FORCE)
SET( APPLICATION_INCLUDES "${INCLUDES}" CACHE STRING "Code-level Application dependency folders" FORCE)
SET( APPLICATION_DEPENDS "${DEPENDENCIES}" CACHE STRING "Code-level Application dependencies" FORCE)
