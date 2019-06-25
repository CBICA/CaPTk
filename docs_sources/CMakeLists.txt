#
# target "Documentation" allows building the Documentationumentation with doxygen/dot on WIN32 and Linux
# Creates .chm windows help file if MS HTML help workshop 
# (available from http://msdn.microsoft.com/workshop/author/htmlhelp)
# is installed with its DLLs in PATH.
#
#
# Please note, that the tools, e.g.:
# doxygen, dot, latex, dvips, makeindex, gswin32, etc.
# must be in path.
#
# Note about Visual Studio Projects: 
# MSVS hast its own path environment which may differ from the shell.
# See "Menu Tools/Options/Projects/VC++ Directories" in VS 7.1
#
# author Jan Woetzel 2004-2006
# www.mip.informatik.uni-kiel.de/~jw


FIND_PACKAGE(Doxygen REQUIRED)

# click+jump in Emacs and Visual Studio (for doxy.config) (jw)
IF    (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv)")
  SET(DOXY_WARN_FORMAT "\"$file($line) : $text \"")
ELSE  (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv)")
  SET(DOXY_WARN_FORMAT "\"$file:$line: $text \"")
ENDIF (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv)")

# we need latex for doxygen because of the formulas
FIND_PACKAGE(LATEX)
IF    (NOT LATEX_COMPILER)
  MESSAGE(STATUS "latex command LATEX_COMPILER not found but usually required. You will probably get warnings and user interaction on doxy run.")
ENDIF (NOT LATEX_COMPILER)
IF    (NOT MAKEINDEX_COMPILER)
  MESSAGE(STATUS "makeindex command MAKEINDEX_COMPILER not found but usually required.")
ENDIF (NOT MAKEINDEX_COMPILER)
IF    (NOT DVIPS_CONVERTER)
  MESSAGE(STATUS "dvips command DVIPS_CONVERTER not found but usually required.")
ENDIF (NOT DVIPS_CONVERTER)

IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/doxy.config.in")
	#MESSAGE(STATUS "configured ${CMAKE_CURRENT_SOURCE_DIR}/doxy.config.in --> ${CMAKE_CURRENT_BINARY_DIR}/doxy.config")
	CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doxy.config.in 
    ${CMAKE_CURRENT_BINARY_DIR}/doxy.config
    @ONLY )
    
	CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/html_footer.html 
    ${CMAKE_CURRENT_BINARY_DIR}/html_footer.html
    @ONLY )
    
	# use (configured) doxy.config from (out of place) BUILD tree:
	SET(DOXY_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/doxy.config")
ELSE (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/doxy.config.in")
	# use static hand-edited doxy.config from SOURCE tree:
	#SET(DOXY_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/doxy.config")
	MESSAGE(SEND_ERROR "Please create ${CMAKE_CURRENT_SOURCE_DIR}/doxy.config.in (or doxy.config as fallback)")
ENDIF( )

FILE(GLOB_RECURSE ALL_DOC_FILES "${PROJECT_SOURCE_DIR}/docs_sources/tutorial/*.*")
FOREACH(docFiles ${ALL_DOC_FILES})
  
  #FILE(COPY ${docFiles}
  #   DESTINATION  ${CMAKE_CURRENT_BINARY_DIR}/html/ )
     
  GET_FILENAME_COMPONENT( TUTORIAL_COMPONENT ${docFiles} NAME_WE )
  
  CONFIGURE_FILE(${docFiles}
    ${CMAKE_CURRENT_BINARY_DIR}/html/${TUTORIAL_COMPONENT}.html
    @ONLY )
     
ENDFOREACH()

FILE(GLOB_RECURSE ALL_IMAGE_FILES "${PROJECT_SOURCE_DIR}/docs_sources/images/*.png")
FOREACH(imgFiles ${ALL_IMAGE_FILES})
  
  FILE(COPY ${imgFiles}
     DESTINATION  ${CMAKE_CURRENT_BINARY_DIR}/html/ )
  
  FILE(COPY ${imgFiles}
     DESTINATION  ${PROJECT_SOURCE_DIR}/docs/ )
  
ENDFOREACH()

IF( ${CAPTK_UPDATE_DOCS} )
  ADD_CUSTOM_TARGET(Documentation ALL ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/html ${PROJECT_SOURCE_DIR}/docs)
ELSE()
  ADD_CUSTOM_TARGET(Documentation ALL ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG} )
ENDIF()
# create a windows help .chm file using hhc.exe
# HTMLHelp DLL must be in path!
# fallback: use hhw.exe interactively
IF(WIN32)
  FIND_PACKAGE(HTMLHelp)
  IF   (HTML_HELP_COMPILER)      
    SET (TMP "${PROJECT__BINARY_DIR}\\Documentation\\html\\index.hhp")
    STRING(REGEX REPLACE "[/]" "\\\\" HHP_FILE ${TMP} )
    # MESSAGE(SEND_ERROR "DBG  HHP_FILE=${HHP_FILE}")
   
    IF (NOT TARGET_DOC_SKIP_INSTALL)
    # install windows help?
    # determine useful name for output file 
    # should be project and version unique to allow installing 
    # multiple projects into one global directory      
    IF   (EXISTS "${CURRENT_BINARY_DIR}/html/index.chm")
      IF   (PROJECT_NAME)
        SET(OUT "${PROJECT_NAME}")
      ELSE (PROJECT_NAME)
        SET(OUT "Documentation") # default
      ENDIF(PROJECT_NAME)
      IF   (${PROJECT_NAME}_VERSION_MAJOR)
        SET(OUT "${OUT}-${${PROJECT_NAME}_VERSION_MAJOR}")
        IF   (${PROJECT_NAME}_VERSION_MINOR)
          SET(OUT  "${OUT}.${${PROJECT_NAME}_VERSION_MINOR}")
          IF   (${PROJECT_NAME}_VERSION_PATCH)
            SET(OUT "${OUT}.${${PROJECT_NAME}_VERSION_PATCH}")      
          ENDIF(${PROJECT_NAME}_VERSION_PATCH)
        ENDIF(${PROJECT_NAME}_VERSION_MINOR)
      ENDIF(${PROJECT_NAME}_VERSION_MAJOR)
      # keep suffix
      SET(OUT  "${OUT}.chm")
      
      #MESSAGE("DBG ${PROJECT_BINARY_DIR}/Documentation/html/index.chm \n${OUT}")
      # create target used by install and package commands 
      INSTALL(FILES "${CURRENT_BINARY_DIR}/html/index.chm"
        DESTINATION Share/docs
        RENAME "${OUT}"
      )
    ENDIF(EXISTS "${CURRENT_BINARY_DIR}/html/index.chm")
    ENDIF(NOT TARGET_DOC_SKIP_INSTALL)

  ENDIF(HTML_HELP_COMPILER)
  # MESSAGE(SEND_ERROR "HTML_HELP_COMPILER=${HTML_HELP_COMPILER}")
ENDIF (WIN32) 

#ADD_CUSTOM_COMMAND( 
#  TARGET Documentation
#  POST_BUILD
#  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex
#  COMMAND make 
#  #OUTPUT refman.pdf
#)

SET( PDF_FILE "${CMAKE_CURRENT_BINARY_DIR}/latex/refman.pdf")

IF(EXISTS ${PDF_FILE})
  CONFIGURE_FILE( 
    ${PDF_FILE}
    ${PROJECT_BINARY_DIR}/Manual.pdf
    COPYONLY
  )
  
  #INSTALL(FILES 
  #  ${PROJECT_BINARY_DIR}/Manual.pdf
  #  DESTINATION ".")

ENDIF()
