SET( PACKAGE_NAME "SBRTLung_Test" )

#PROJECT( ${PACKAGE_NAME} )

INCLUDE_DIRECTORIES( "${CMAKE_SOURCE_DIR}/../common_includes" )


INCLUDE(${ITK_USE_FILE})

# remove some MSVC flags
IF(MSVC)
	ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
ENDIF(MSVC)

ADD_EXECUTABLE( ${PACKAGE_NAME}
  ${Common_dependents}
  ${segment_src}
  ${analyze_src}
  ${CMAKE_CURRENT_SOURCE_DIR}/sbrtlung_test.cxx
  
)

TARGET_LINK_LIBRARIES(
  ${PACKAGE_NAME}
  ${ITK_LIBRARIES}

)

# Test for working memory buffer
#ADD_TEST( NAME Buffer_Test COMMAND ${PACKAGE_NAME} -b "test" )

# Test for CmdParser class
#ADD_TEST( NAME CmdParser_Test COMMAND ${PACKAGE_NAME} -cmdParser "random" )

# Test for Logger class
#ADD_TEST( NAME Logger_Test COMMAND ${PACKAGE_NAME} -logging "testing_phrase" )

# Test for CSV Parser
#ADD_TEST( NAME CSVParser_Test COMMAND ${PACKAGE_NAME} -csvParser ${CMAKE_CURRENT_SOURCE_DIR}/list.csv )

# Test for Geodesic class
ADD_TEST( NAME SBRTlung_Test COMMAND ${PACKAGE_NAME} -g)