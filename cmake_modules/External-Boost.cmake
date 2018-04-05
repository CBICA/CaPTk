## DO NOT USE!!
## This is an incomplete script. Boost currently does NOT support CMake as a building platform. 
## For more details, visit: https://svn.boost.org/trac/boost/wiki/CMakeModularizationStatus 

MESSAGE( "External project - Boost" )

SET( Boost_DEPENDENCIES )

ExternalProject_Add( Boost
	DEPENDS ${Boost_DEPENDENCIES}
	SVN_REPOSITORY "https://svn.boost.org/svn/boost/tags/release/Boost_1_42_0/"
	SOURCE_DIR Boost
	BINARY_DIR Boost-build
	UPDATE_COMMAND ""
	PATCH_COMMAND ""
	CMAKE_GENERATOR ${gen}
	CMAKE_ARGS
		${ep_common_args}
		-DBUILD_EXAMPLES:BOOL=OFF
		-DBUILD_SHARED_LIBS:BOOL=ON
		-DBUILD_TESTING:BOOL=OFF
		-DCMAKE_BUILD_TYPE:STRING=Release
		-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
		-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/INSTALL
	)

set( Boost_DIR ${CMAKE_BINARY_DIR}/Boost-build )