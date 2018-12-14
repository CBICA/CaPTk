# Macro to add, link, build, and install an application and it's cwl file
MACRO(CAPTK_ADD_EXECUTABLE APPLICATION)
  ADD_EXECUTABLE(${APPLICATION}
    ${PROJECT_SOURCE_DIR}/src/applications/${APPLICATION}.cxx
    #${CBICA-TK_SOURCES}
    #${YAMLCPP_Headers}
    #${YAMLCPP_Sources}
    #${APPLICATION_DEPENDS}
  )

  TARGET_LINK_LIBRARIES(${APPLICATION}
    #${DEPENDENT_LIBS}
    ${LIBNAME_CBICATK}
    ${LIBNAME_Applications}
    ${LIBNAME_FeatureExtractor}
  )

  SET_TARGET_PROPERTIES( ${APPLICATION} PROPERTIES FOLDER "${StandAloneCLIAppsFolder}" )

  ADD_DEPENDENCIES( ${APPLICATION} ${LIBNAME_Applications} ${LIBNAME_FeatureExtractor} ${LIBNAME_CBICATK} )

  IF (APPLE) 
    # list (APPEND STANDALONE_APPS_LIST ${APPLICATION})
    INSTALL( TARGETS ${APPLICATION}
    BUNDLE DESTINATION .
    RUNTIME DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
    LIBRARY DESTINATION ${EXE_NAME}.app/Contents/Resources/lib
    CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  ) 

  ELSE()
    INSTALL( TARGETS ${APPLICATION}
    BUNDLE DESTINATION .
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    CONFIGURATIONS "${CMAKE_CONFIGURATION_TYPES}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )

  ENDIF()

  # Add test for run tests
  ADD_TEST( NAME ${APPLICATION}_rt COMMAND ${APPLICATION} -rt )

  CWL_INSTALL(${APPLICATION})
  
ENDMACRO()

# Macro to generate and install a cwl file after target application is built
MACRO(CWL_INSTALL APPLICATION)

  # Post build cwl generation
  add_custom_command(TARGET ${APPLICATION}
    POST_BUILD
    COMMAND ${APPLICATION} -cwl
    COMMENT "Generating cwl for ${APPLICATION}..."
    VERBATIM
  )

  IF (APPLE) 
    # list (APPEND STANDALONE_APPS_LIST ${APPLICATION})
    INSTALL( FILES ${PROJECT_BINARY_DIR}/${APPLICATION}.cwl
    DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
  ) 

  ELSE()
    INSTALL( FILES ${PROJECT_BINARY_DIR}/${APPLICATION}.cwl
    DESTINATION bin
  )

  ENDIF()

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
